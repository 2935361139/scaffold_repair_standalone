import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rdkit import Chem
from scaffold_repair.import_state import block,restore
from universal_input.normalize import convert,build_graph,heavy_view,read_record

class TestUniversalInput(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.ref=Chem.MolFromSmiles('CC')
        conf=Chem.Conformer(2);conf.Set3D(True);conf.SetAtomPosition(0,(0,0,0));conf.SetAtomPosition(1,(1.5,0,0));self.ref.AddConformer(conf)
        (self.p/'reference.sdf').write_text(block(self.ref))
        (self.p/'protein.pdb').write_text('ATOM      1  CA  ALA A   1      20.000  20.000  20.000  1.00 20.00           C  \nEND\n')
    def tearDown(self):self.tmp.cleanup()
    def recipe(self, output='out', protein=True):
        return dict(format='sdf',source='input.sdf',mapping_file='mapping.json',reference_sdf='reference.sdf',output=output,frame='reference',protein_pdb='protein.pdb' if protein else '',protein_frame_confirmed=protein,experiment_output='experiment')
    def write_sdf_input(self):
        (self.p/'input.sdf').write_text(block(self.ref))
        (self.p/'mapping.json').write_text(json.dumps({'samples':[{'id':'s1','sdf_index':0,'fixed_generated_indices':[0],'fixed_reference_indices':[0]}]}))
    def test_sdf_ready_and_v2_reader(self):
        self.write_sdf_input();r=convert(self.recipe(),self.p);self.assertEqual(r['status'],'ready')
        from scaffold_repair.adapters.sdf import load_manifest
        rows=load_manifest(self.p/'out/manifest.json');self.assertEqual(len(rows),1);self.assertEqual(rows[0]['mol'].GetNumAtoms(),2)
    def test_missing_protein_pending(self):
        self.write_sdf_input();r=convert(self.recipe(protein=False),self.p)
        self.assertEqual(r['converted_samples'],1);self.assertFalse((self.p/'out/run_config.py').exists())
    def test_author_labels_do_not_change_structure(self):
        self.write_sdf_input();c=self.recipe('a');c['model']='FLOWR';convert(c,self.p)
        m=json.loads((self.p/'mapping.json').read_text());m['samples'][0]['success']=False;m['samples'][0]['cause']='anything';(self.p/'mapping.json').write_text(json.dumps(m))
        c=self.recipe('b');c['model']='unseen_generator';convert(c,self.p)
        self.assertEqual((self.p/'a/samples/s1/molecule.sdf').read_bytes(),(self.p/'b/samples/s1/molecule.sdf').read_bytes())
    def test_refuses_partial_mapping(self):
        self.write_sdf_input();(self.p/'input.sdf').write_text(block(self.ref)*2)
        with self.assertRaises(ValueError):convert(self.recipe(),self.p)
    def test_overwrite_blocked(self):
        self.write_sdf_input();convert(self.recipe(),self.p)
        with self.assertRaises(ValueError):convert(self.recipe(),self.p)
    def test_custom_npz_fields(self):
        np.savez(self.p/'graphs.npz',species=np.array([6,6]),q=np.array([0,0]),adj=np.array([[0,7],[7,0]]),xyz=np.array([[0,0,0],[.15,0,0]]))
        (self.p/'mapping.json').write_text(json.dumps({'samples':[{'id':'s1','fixed_generated_indices':[0],'fixed_reference_indices':[0]}]}))
        c=self.recipe();c.update(format='npz_graphs',source='graphs.npz',units='nanometer',graph_profile={'keys':{'atoms':'species','charges':'q','bonds':'adj','positions':'xyz'},'atom_encoding':'atomic_number','bond_decoder':{'7':'SINGLE'}})
        r=convert(c,self.p);self.assertEqual(r['status'],'ready');self.assertAlmostEqual(np.load(self.p/'out/samples/s1/coordinates.npy')[1,0],1.5)
    def test_asymmetric_npz_rejected(self):
        a={'atoms':np.array([6,6]),'charges':np.zeros(2),'bonds':np.array([[0,1],[0,0]]),'positions':np.zeros((2,3))}
        with self.assertRaises(ValueError):build_graph(a,[],{'1':'SINGLE'},'atomic_number')
    def test_overvalence_survives_hydrogen_removal(self):
        a={'atoms':np.array([6,1,1,1,1,1]),'charges':np.zeros(6),'bonds':np.zeros((6,6),int),'positions':np.zeros((6,3))}
        a['bonds'][0,1:]=1;a['bonds'][1:,0]=1
        raw=build_graph(a,[],{'1':'SINGLE'},'atomic_number');h,m=heavy_view(raw);h.UpdatePropertyCache(strict=False)
        parsed=restore(next(Chem.ForwardSDMolSupplier(io.BytesIO(block(h).encode()),sanitize=False,removeHs=False)))
        self.assertEqual(parsed.GetAtomWithIdx(0).GetNumExplicitHs(),5)
        with self.assertRaises(Exception):Chem.SanitizeMol(parsed)
    def test_nonstandard_hydrogen_rejected(self):
        raw=Chem.MolFromSmiles('[H-]')
        with self.assertRaises(ValueError):heavy_view(raw)
    def test_aromatic_nh_survives_serialization(self):
        m=Chem.MolFromSmiles('c1cc[nH]c1');m.UpdatePropertyCache(strict=False)
        parsed=restore(next(Chem.ForwardSDMolSupplier(io.BytesIO(block(m).encode()),sanitize=False,removeHs=False)))
        Chem.SanitizeMol(parsed)
        self.assertEqual(Chem.MolToSmiles(m),Chem.MolToSmiles(parsed))
    def test_bad_mapping_does_not_disappear(self):
        self.write_sdf_input();p=self.p/'mapping.json';m=json.loads(p.read_text());m['samples'][0]['fixed_generated_indices']=[99];p.write_text(json.dumps(m))
        r=convert(self.recipe(),self.p);self.assertEqual(r['total_samples'],1);self.assertEqual(len(r['conversion_errors']),1);self.assertFalse((self.p/'out/run_config.py').exists())
    def test_readwrite_preserves_atom_state(self):
        from scaffold_repair.evaluation import read_mol,save_mol
        m=Chem.MolFromSmiles('c1cc[nH]c1');p=self.p/'saved.sdf';p.write_text(block(m));mol=read_mol(p);save_mol(self.p/'resaved.sdf',mol);again=read_mol(self.p/'resaved.sdf');Chem.SanitizeMol(again)
        self.assertEqual(Chem.MolToSmiles(m),Chem.MolToSmiles(again))

    def test_bond_end_order_is_not_a_different_graph(self):
        from scaffold_repair.import_state import encode,PROPERTY
        m=encode(Chem.Mol(self.ref));state=json.loads(m.GetProp(PROPERTY))
        state['bonds']=[[b,a,k] for a,b,k in state['bonds']]
        m.SetProp(PROPERTY,json.dumps(state));restore(m)
        state['bonds'][0][2]='DOUBLE';m.SetProp(PROPERTY,json.dumps(state))
        with self.assertRaises(ValueError):restore(m)

    def test_chiral_disconnected_sdf_keeps_structure(self):
        m=Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)C(=O)O.[Na+]'))
        from rdkit.Chem import AllChem
        AllChem.EmbedMolecule(m,randomSeed=7)
        heavy,_=heavy_view(m)
        parsed=restore(next(Chem.ForwardSDMolSupplier(io.BytesIO(block(heavy).encode()),sanitize=False,removeHs=False)))
        Chem.SanitizeMol(parsed)
        self.assertEqual(Chem.MolToSmiles(heavy),Chem.MolToSmiles(parsed))
        self.assertEqual(len(Chem.GetMolFrags(parsed)),2)

if __name__=='__main__':unittest.main()

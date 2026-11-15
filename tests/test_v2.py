"""Model-neutral graph, receptor and accounting invariants; no case-name rules."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from rdkit import Chem
from scaffold_repair.topology import plan, LIMITS, fixed_matches, graph_key, graph_edits
from scaffold_repair.receptor import prepare
from scaffold_repair.geometry_intervention import route_failure
from scaffold_repair.evaluation import task_from_mapping
from posebusters.tools.protein import get_atom_type_mask


def molecule(smiles, coords):
    m=Chem.MolFromSmiles(smiles)
    c=Chem.Conformer(m.GetNumAtoms())
    for i,x in enumerate(coords): c.SetAtomPosition(i,x)
    m.AddConformer(c)
    return m


def item(mol,ref,fixed,ref_fixed=None):
    return dict(mol=mol,ref=ref,slots=list(range(mol.GetNumAtoms())),
                task=task_from_mapping(ref,fixed,ref_fixed or fixed),id='arbitrary_model_sample')


def pdb_atom(mol,idx,name,number,hetero=True,element_name=''):
    info=Chem.AtomPDBResidueInfo()
    info.SetName((element_name or mol.GetAtomWithIdx(idx).GetSymbol()).rjust(4))
    info.SetResidueName(name);info.SetResidueNumber(number);info.SetChainId('A')
    info.SetIsHeteroAtom(hetero);info.SetSerialNumber(idx+1)
    mol.GetAtomWithIdx(idx).SetMonomerInfo(info)


class TopologyTests(unittest.TestCase):
    def test_missing_bridge_is_proposed_and_atoms_are_preserved(self):
        ref=molecule('CCO',[(0,0,0),(1.5,0,0),(2.9,.1,0)])
        broken=Chem.RWMol(ref);broken.RemoveBond(0,1);broken=broken.GetMol();Chem.SanitizeMol(broken)
        proposals,audit=plan(item(broken,ref,[0,2]))
        self.assertTrue(any(c['family']=='L' and len(Chem.GetMolFrags(c['mol']))==1 for c in proposals))
        for c in proposals:
            self.assertEqual(c['mol'].GetNumAtoms(),3)
            np.testing.assert_array_equal(c['mol'].GetConformer().GetPositions(),broken.GetConformer().GetPositions())
        self.assertLessEqual(audit['attempted_sanitizations'],LIMITS['max_sanitizations'])

    def test_fixed_double_bond_is_restored_without_copying_core(self):
        ref=molecule('CC(=O)N',[(0,0,0),(1.5,0,0),(1.5,1.2,0),(2.8,0,0)])
        broken=Chem.RWMol(ref);broken.GetBondBetweenAtoms(1,2).SetBondType(Chem.BondType.SINGLE)
        broken=broken.GetMol();Chem.SanitizeMol(broken)
        proposals,_=plan(item(broken,ref,[1,2]))
        self.assertTrue(proposals)
        self.assertTrue(all(fixed_matches(c['mol'],ref,{1:1,2:2}) for c in proposals if c['family']=='F'))
        self.assertEqual(proposals[0]['edits']['bonds'],[dict(atoms=[1,2],before='SINGLE',after='DOUBLE')])

    def test_incompatible_external_halogen_bond_is_not_left_overvalent(self):
        ref=molecule('FCCC',[(0,0,0),(1.4,0,0),(2.9,0,0),(4.4,0,0)])
        broken=Chem.RWMol(ref);broken.RemoveBond(0,1);broken.AddBond(0,3,Chem.BondType.SINGLE)
        broken=broken.GetMol();Chem.SanitizeMol(broken)
        proposals,_=plan(item(broken,ref,[0,1]))
        good=[c for c in proposals if c['family']=='F']
        self.assertTrue(good)
        for c in good:
            self.assertIsNotNone(c['mol'].GetBondBetweenAtoms(0,1))
            self.assertIsNone(c['mol'].GetBondBetweenAtoms(0,3))
            Chem.SanitizeMol(c['mol'])

    def test_aromatic_serialization_is_not_a_repair(self):
        ref=molecule('c1ccccc1C',[(0,0,0),(1,0,0),(2,1,0),(2,2,0),(1,3,0),(0,2,0),(-1,0,0)])
        kek=Chem.Mol(ref);Chem.Kekulize(kek,clearAromaticFlags=True)
        proposals,_=plan(item(kek,ref,list(range(6))))
        self.assertEqual(proposals,[])
        self.assertEqual(graph_edits(kek,ref)['bonds'],[])

    def test_far_components_are_not_falsely_connected(self):
        ref=molecule('CC',[(0,0,0),(20,0,0)])
        broken=Chem.RWMol(ref);broken.RemoveBond(0,1);broken=broken.GetMol();Chem.SanitizeMol(broken)
        proposals,audit=plan(item(broken,ref,[0]))
        self.assertEqual(proposals,[])
        self.assertEqual(audit['status'],'no_admissible_candidate_within_budget')

    def test_fixed_fixed_new_bridge_is_forbidden(self):
        ref=molecule('C.C',[(0,0,0),(1.5,0,0)])
        proposals,_=plan(item(ref,ref,[0,1]))
        self.assertEqual(proposals,[])

    def test_metadata_does_not_select_model_specific_rules(self):
        ref=molecule('CCO',[(0,0,0),(1.5,0,0),(2.9,0,0)])
        broken=Chem.RWMol(ref);broken.RemoveBond(0,1);broken=broken.GetMol();Chem.SanitizeMol(broken)
        a=item(broken,ref,[0,2]);b=dict(a,id='a_completely_new_dataset',model='unseen_generator')
        pa,aa=plan(a);pb,ab=plan(b)
        self.assertEqual([graph_key(c['mol']) for c in pa],[graph_key(c['mol']) for c in pb])
        self.assertEqual(aa,ab)

    def test_unsupported_inputs_are_explicit(self):
        row=dict(molecule_readable=False,rdkit_valid=False)
        self.assertEqual(route_failure(row)['status'],'unsupported')
        row=dict(molecule_readable=True,rdkit_valid=True,fixed_element_changed_sites=[2])
        self.assertEqual(route_failure(row)['status'],'unsupported')


class ReceptorTests(unittest.TestCase):
    def receptor(self):
        # Organic chloride belongs to its organic residue, not to metal/ion class.
        env=molecule('C.CCl.O.[Cl-]',[(10,0,0),(0,0,0),(1.7,0,0),(4,0,0),(6,0,0)])
        pdb_atom(env,0,'ALA',1,False)
        pdb_atom(env,1,'LIG',2,True,'C1');pdb_atom(env,2,'LIG',2,True,'CL1')
        pdb_atom(env,3,'HOH',3);pdb_atom(env,4,'CL',4)
        return env

    def test_environment_categories_are_disjoint(self):
        env=self.receptor(); kinds={'protein','waters','organic_cofactors','inorganic_cofactors'}
        expected={'protein':[0],'organic_cofactors':[1,2],'waters':[3],'inorganic_cofactors':[4]}
        for kind,indices in expected.items():
            mask=get_atom_type_mask(env,kinds-{kind}|{'hydrogens'})
            self.assertEqual([i for i,b in enumerate(mask) if b],indices)

    def test_reference_exclusion_preserves_water_and_original(self):
        env=self.receptor();ref=molecule('CCl',[(0,0,0),(1.7,0,0)])
        with tempfile.TemporaryDirectory() as t:
            source=Path(t)/'complex.pdb';dest=Path(t)/'receptor.pdb'
            Chem.MolToPDBFile(env,str(source));before=source.read_bytes()
            audit=prepare(source,ref,dest)
            self.assertEqual(source.read_bytes(),before)
            self.assertEqual(audit['automatic_exclusions'],[['A',2,'','LIG']])
            result=Chem.MolFromPDBFile(str(dest),sanitize=False,removeHs=False,proximityBonding=False)
            names=[a.GetPDBResidueInfo().GetResidueName().strip() for a in result.GetAtoms()]
            self.assertIn('HOH',names);self.assertIn('ALA',names);self.assertNotIn('LIG',names)

    def test_unmatched_cofactor_is_not_removed(self):
        env=self.receptor();ref=molecule('CCl',[(20,0,0),(21.7,0,0)])
        with tempfile.TemporaryDirectory() as t:
            source=Path(t)/'complex.pdb';dest=Path(t)/'receptor.pdb'
            Chem.MolToPDBFile(env,str(source))
            audit=prepare(source,ref,dest)
            self.assertEqual(audit['automatic_exclusions'],[])
            self.assertEqual(source.read_bytes(),dest.read_bytes())


if __name__=='__main__':unittest.main()

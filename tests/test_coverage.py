"""Independent failure classes and bounded intervention contracts."""
import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from test_v2 import molecule,item
from scaffold_repair.chemistry import problems
from scaffold_repair.topology import plan,LIMITS,graph_key
from scaffold_repair.conformers import propose,SETTINGS
from scaffold_repair.coverage import fixed_stereo,batch_metrics,geometry_details,CATALOG
from posebusters import PoseBusters


class CoverageTests(unittest.TestCase):
    def test_generated_region_overvalence_has_a_bounded_generic_hypothesis(self):
        bad=Chem.MolFromSmiles('C(C)(C)(C)=C',sanitize=False)
        conf=Chem.Conformer(5)
        for i in range(5): conf.SetAtomPosition(i,(i*1.4,0.,0.))
        bad.AddConformer(conf)
        ref=molecule('CC(C)(C)C',[(i*1.4,0.,0.) for i in range(5)])
        self.assertTrue(problems(bad))
        p,a=plan(item(bad,ref,[2],[0]))
        self.assertTrue(any(c['family']=='E' and len(Chem.GetMolFrags(c['mol']))==1 for c in p))
        self.assertLessEqual(a['attempted_sanitizations'],LIMITS['max_sanitizations'])
        self.assertLessEqual(len(p),LIMITS['max_candidates'])

    def test_conformer_rebuild_preserves_graph_and_fixed_coordinates(self):
        m=Chem.AddHs(Chem.MolFromSmiles('CCCO'))
        self.assertEqual(AllChem.EmbedMolecule(m,randomSeed=31),0)
        m=Chem.RemoveHs(m)
        x=m.GetConformer().GetPositions()
        inp=item(m,m,[0,1]);inp['xyz']=x;inp['target']=x.copy()
        candidates,audit=propose(inp,[])
        self.assertTrue(candidates)
        self.assertLessEqual(len(candidates),SETTINGS['max_candidates'])
        for c in candidates:
            self.assertEqual(graph_key(c['mol']),graph_key(m))
            np.testing.assert_array_equal(c['coordinates'][[0,1]],x[[0,1]])
        renamed=dict(inp,id='unseen_model_100')
        again,_=propose(renamed,[])
        for a,b in zip(candidates,again): np.testing.assert_array_equal(a['coordinates'],b['coordinates'])

    def test_disconnected_graph_is_not_disguised_as_a_new_conformer(self):
        m=molecule('CC.C',[(0.,0.,0.),(1.5,0.,0.),(20.,0.,0.)])
        inp=item(m,m,[0]);inp['xyz']=m.GetConformer().GetPositions()
        p,a=propose(inp,[])
        self.assertEqual(p,[])
        self.assertEqual(a['status'],'no_eligible_connected_graph')

    def test_local_stereo_inversion_is_reported_without_a_cip_shortcut(self):
        ref=molecule('F[C@](Cl)(Br)I',[(1.,0.,0.),(0.,0.,0.),(0.,1.,0.),(0.,0.,1.),(-1.,-1.,-1.)])
        changed=Chem.Mol(ref)
        for i,x in enumerate(ref.GetConformer().GetPositions()): changed.GetConformer().SetAtomPosition(i,(-x[0],x[1],x[2]))
        r=fixed_stereo(item(changed,ref,list(range(5))),changed)
        self.assertTrue(any(c['status']=='inverted_orientation' for c in r['centers']))

    def test_invalid_bond_geometry_is_localized_to_atoms(self):
        m=molecule('CCO',[(0.,0.,0.),(5.,0.,0.),(6.4,.1,0.)])
        loc=geometry_details(m,PoseBusters('dock').config)
        self.assertTrue(any(r['type']=='bonds' and r['atoms']==[0,1] for r in loc))

    def test_batch_metrics_do_not_call_missing_novelty_a_pass(self):
        m=batch_metrics(['CCO','CCO','CC',None],4)
        self.assertEqual(m['rdkit_valid_records'],3)
        self.assertEqual(m['unique_molecules'],2)
        self.assertEqual(m['training_novelty']['status'],'not_assessed')
        self.assertEqual(m['diversity_pairs'],3)

    def test_catalog_separates_structural_repairs_from_drug_discovery_oracles(self):
        cats={c['id']:c for c in CATALOG}
        self.assertIn('E',cats['chemical_graph']['interventions'])
        self.assertIn('G',cats['intramolecular_geometry']['interventions'])
        self.assertEqual(cats['binding']['interventions'],[])
        self.assertEqual(cats['distribution']['interventions'],[])


if __name__=='__main__': unittest.main()

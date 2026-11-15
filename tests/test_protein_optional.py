"""Scope correctness, strict input errors, unchanged protein checks and scoped replay."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from posebusters import PoseBusters
from scaffold_repair.context import (protein_input, molecular_config, omitted_checks,
                                     MOLECULAR_SCOPE, PROTEIN_SCOPE)
from scaffold_repair.adapters.sdf import load_one, load_manifest
from scaffold_repair.evaluation import evaluate, task_from_mapping
from scaffold_repair.legacy_functions import pb_metrics
from scaffold_repair.pipeline import check, intervene, restore, verify_checked
from scaffold_repair.engine import plan_bank, ARMS, MOLECULAR_ARMS
from scaffold_repair.import_state import block
from universal_input.normalize import convert


class ProteinOptionalTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)
        mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
        self.assertEqual(AllChem.EmbedMolecule(mol, randomSeed=17), 0)
        AllChem.UFFOptimizeMolecule(mol)
        self.ref = Chem.RemoveHs(mol)
        (self.p/'reference.sdf').write_text(block(self.ref))
        (self.p/'generated.sdf').write_text(block(self.ref))
        (self.p/'protein.pdb').write_text('ATOM      1  CA  ALA A   1       4.000   4.000   4.000  1.00 20.00           C  \nEND\n')
        self.spec = dict(id='sample', input_id='task', sdf='generated.sdf', reference_sdf='reference.sdf',
                         protein_pdb='', protein_mode='absent', fixed_generated_indices=[0], fixed_reference_indices=[0])

    def tearDown(self):
        self.tmp.cleanup()

    def evaluation(self, mol=None, protein=None):
        mol = self.ref if mol is None else mol
        return evaluate(mol, mol.GetConformer().GetPositions(), self.ref,
                        task_from_mapping(self.ref, [0], [0]), list(range(mol.GetNumAtoms())), protein)[0]

    def recipe(self):
        (self.p/'mapping.json').write_text(json.dumps({'samples':[
            {'id':'sample','sdf_index':0,'fixed_generated_indices':[0],'fixed_reference_indices':[0]}]}))
        return dict(format='sdf',source='generated.sdf',reference_sdf='reference.sdf',mapping_file='mapping.json',
                    output='converted',experiment_output='experiment',frame='reference',protein_mode='absent',protein_pdb='')

    def test_missing_protein_does_not_opt_in(self):
        spec = dict(self.spec); del spec['protein_mode']
        with self.assertRaisesRegex(ValueError, 'Protein is required'):
            load_one(spec, self.p)

    def test_broken_pdb_does_not_fall_back(self):
        with self.assertRaisesRegex(ValueError, 'missing'):
            load_one(dict(self.spec, protein_mode='required', protein_pdb='missing.pdb'), self.p)
        with self.assertRaisesRegex(ValueError, 'conflicts'):
            load_one(dict(self.spec, protein_pdb='missing.pdb'), self.p)

    def test_unknown_mode_rejected(self):
        with self.assertRaisesRegex(ValueError, 'protein_mode'):
            protein_input('automatic', '', self.p)

    def test_malformed_pdb_is_runtime_input_error(self):
        (self.p/'bad.pdb').write_text('invalid receptor\n')
        item = load_one(dict(self.spec, protein_mode='required', protein_pdb='bad.pdb'), self.p)
        with self.assertRaisesRegex(ValueError, 'Evaluation stopped'):
            check([item], self.p/'bad_run')
        summary = json.loads((self.p/'bad_run/summary.json').read_text())
        self.assertFalse(summary['completed']); self.assertEqual(summary['evaluated_samples'], 0)

    def test_conversion_absent_mode_produces_executable_manifest(self):
        result = convert(self.recipe(), self.p)
        self.assertEqual(result['status'], 'ready')
        self.assertTrue((self.p/'converted/run_config.py').is_file())
        item = load_manifest(self.p/'converted/manifest.json')[0]
        self.assertIsNone(item['protein']); self.assertEqual(item['protein_mode'], 'absent')

    def test_converter_conflict_is_not_a_partial_batch(self):
        cfg = dict(self.recipe(), protein_pdb='missing.pdb')
        result = convert(cfg, self.p)
        self.assertEqual(len(result['conversion_errors']), 1)
        self.assertFalse((self.p/'converted/run_config.py').exists())

    def test_molecular_thresholds_are_identical_subsets(self):
        original = PoseBusters('dock').config
        subset = molecular_config(original)
        for module in subset['modules']:
            if module['function'] != 'loading':
                self.assertIn(module, original['modules'])
        self.assertEqual(len(omitted_checks(original)), 10)
        self.assertIn('mol_cond_loaded', omitted_checks(original))

    def test_absent_environment_is_null_never_passed(self):
        row = self.evaluation()
        self.assertEqual(row['evaluation_scope'], MOLECULAR_SCOPE)
        self.assertTrue(row['posebusters_valid'])
        self.assertTrue(row['not_evaluated_checks'])
        for key, detail in row['not_evaluated_checks'].items():
            self.assertNotIn(key, row['posebusters_checks'])
            self.assertIsNone(detail['passed'])
            self.assertEqual(detail['status'], 'not_evaluated')

    def test_bad_bond_is_still_detected_without_protein(self):
        mol = Chem.Mol(self.ref)
        mol.GetConformer().SetAtomPosition(2, (20.,0.,0.))
        row = self.evaluation(mol)
        self.assertFalse(row['posebusters_valid'])
        self.assertIn('bond_lengths', row['posebusters_failed_checks'])

    def test_invalid_chemistry_keeps_scope_and_is_not_passed(self):
        mol = Chem.Mol(self.ref); mol.GetAtomWithIdx(0).SetNumExplicitHs(5)
        row = self.evaluation(mol)
        self.assertFalse(row['rdkit_valid'])
        self.assertEqual(row['evaluation_scope'], MOLECULAR_SCOPE)
        self.assertTrue(row['not_evaluated_checks'])

    def test_with_protein_native_evaluation_is_unchanged(self):
        protein = self.p/'protein.pdb'
        old = pb_metrics(self.ref, protein, PoseBusters, Chem)
        new = self.evaluation(protein=protein)
        for key, value in old.items(): self.assertEqual(new[key], value, key)
        self.assertEqual(new['evaluation_scope'], PROTEIN_SCOPE)
        self.assertEqual(new['not_evaluated_checks'], {})

    def test_scoped_snapshot_replay_and_no_fake_protein(self):
        item = load_one(self.spec, self.p)
        result = check([item], self.p/'generated')
        self.assertTrue(result['completed']); self.assertEqual(result['passed'], 1)
        restored = restore(self.p/'generated/samples/sample')
        self.assertIsNone(restored['protein'])
        self.assertFalse((self.p/'generated/samples/sample/protein.pdb').exists())
        outcome = intervene(self.p/'generated/failed', self.p/'intervention')
        self.assertTrue(outcome['completed'])
        self.assertEqual(outcome['by_evaluation_scope'][MOLECULAR_SCOPE]['passed'], 1)

    def test_receptor_scope_tampering_cannot_enter_intervention(self):
        check([load_one(self.spec,self.p)], self.p/'generated')
        path = self.p/'generated/samples/sample/sample.json'
        data = json.loads(path.read_text()); data['protein_mode']='required'
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'evidence changed'):
            verify_checked(self.p/'generated')

    def test_same_input_id_cannot_mix_contexts(self):
        absent = load_one(self.spec,self.p)
        present = load_one(dict(self.spec,id='second',protein_mode='required',protein_pdb='protein.pdb'),self.p)
        with self.assertRaisesRegex(ValueError,'inconsistent conditions'):
            check([absent,present],self.p/'mixed')

    def test_mixed_batch_reports_separate_scopes(self):
        absent = load_one(self.spec,self.p)
        present = load_one(dict(self.spec,id='second',input_id='protein_task',protein_mode='required',protein_pdb='protein.pdb'),self.p)
        result = check([absent,present],self.p/'mixed')
        self.assertEqual(result['by_evaluation_scope'][MOLECULAR_SCOPE]['total'],1)
        self.assertEqual(result['by_evaluation_scope'][PROTEIN_SCOPE]['total'],1)

    def test_protein_free_plan_excludes_duplicate_environment_controls(self):
        item = load_one(self.spec,self.p)
        plan = plan_bank(item,self.evaluation())
        self.assertEqual(plan['arm_names'],MOLECULAR_ARMS)
        self.assertNotIn('UFF_only',plan['arm_names'])
        item['protein']=self.p/'protein.pdb'
        self.assertEqual(plan_bank(item,self.evaluation())['arm_names'],ARMS)


if __name__ == '__main__': unittest.main()

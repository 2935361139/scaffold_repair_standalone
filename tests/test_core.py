"""Fast invariants; run with: python -m unittest discover -s tests -v"""
import ast
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from rdkit import Chem
from scaffold_repair.geometry_intervention import (repulsion, bounded_minimize, route_failure,
                                                  select_result, successful)
from scaffold_repair.legacy_functions import coordinate_arm, suffix
from scaffold_repair.evaluation import task_from_mapping, fixed_substructure_match, expected_fixed
from scaffold_repair.adapters.flowr_records import recover_float32_origin
from scaffold_repair.adapters.sdf import load_manifest
from scaffold_repair.io import safe_id
from scaffold_repair.pipeline import check


def passing():
    return dict(rdkit_valid=True, fc_valid=True, root_fixed_substructure_match=True,
                fixed_chemistry_preserved=True, posebusters_valid=True, posebusters_failed_checks=[])


class ScientificInvariants(unittest.TestCase):
    def test_primary_priority_and_relaxed_constraint(self):
        base = passing(); base['posebusters_valid'] = False
        self.assertEqual(select_result(base, {'A': passing(), 'C': passing()})['all_successful_interventions'], ['A', 'C'])
        self.assertEqual(select_result(base, {'A': passing(), 'C': passing()})['selected_arm'], 'A')
        self.assertTrue(select_result(base, {'C': passing()})['fixed_coordinate_relaxation'])

    def test_controls_cannot_supply_a_repair(self):
        base = passing(); base['posebusters_valid'] = False
        self.assertIsNone(select_result(base, {'UFF_only': passing(), 'A_random_0': passing()})['selected_arm'])

    def test_original_success_is_preserved(self):
        self.assertEqual(select_result(passing(), {'C': passing()})['selected_arm'], 'baseline')

    def test_discrete_failure_is_unsupported(self):
        base = passing(); base['fc_valid'] = False
        self.assertEqual(route_failure(base)['status'], 'unsupported')
        self.assertFalse(successful(base))

    def test_geometry_routing_is_not_water_specific(self):
        for label in ('minimum_distance_to_protein', 'minimum_distance_to_waters', 'bond_angles', 'internal_energy'):
            base = passing(); base.update(posebusters_valid=False, posebusters_failed_checks=[label])
            self.assertEqual(route_failure(base)['status'], 'supported_geometry')

    def test_unknown_check_is_not_silently_accepted(self):
        base = passing(); base['posebusters_valid'] = None
        self.assertFalse(successful(base))
        self.assertEqual(route_failure(base)['status'], 'unsupported')

    def test_random_control_preserves_magnitude_scope_and_dtype(self):
        start = np.array([[1, 2, 3], [4, 5, 6], [9, 9, 9]], dtype=np.float32)
        goal = start.copy(); goal[0] += [1, 0, 0]; goal[1] += [0, 2, 0]
        result, _ = coordinate_arm(start, goal, [0, 1], 'random_direction', 27)
        self.assertEqual(result.dtype, start.dtype)
        np.testing.assert_array_equal(result[2], start[2])
        np.testing.assert_allclose(np.linalg.norm(result[:2]-start[:2], axis=1), [1, 2], atol=1e-6)
        np.testing.assert_allclose(np.sum((result[:2]-start[:2])*(goal[:2]-start[:2]), axis=1), 0, atol=1e-6)

    def test_repulsion_gradient(self):
        x = np.array([[1.2, .2, .3]]); env = np.array([[0., 0., 0.]])
        cutoff = np.array([[2.5]])
        _, g = repulsion(x, env, cutoff, 500.)
        for k in range(3):
            d = np.zeros_like(x); d[0, k] = 1e-6
            numerical = (repulsion(x+d, env, cutoff, 500.)[0] - repulsion(x-d, env, cutoff, 500.)[0]) / 2e-6
            self.assertAlmostEqual(numerical, g[0, k], places=5)

    def test_optimizer_hard_budget(self):
        result, receipt = bounded_minimize(lambda x: (float(np.dot(x, x)), 2*x),
                                           np.array([1.]), [(-2, 2)], 1, 1)
        self.assertLessEqual(receipt['objective_calls'], 1)
        self.assertLessEqual(receipt['iterations'], 1)
        self.assertTrue(np.isfinite(result).all())

    def test_atom_mapping_rejects_duplicates(self):
        with self.assertRaises(ValueError):
            task_from_mapping(Chem.MolFromSmiles('CCO'), [0, 0], [0, 1])

    def test_aromatic_fixed_query_remains_strict(self):
        aromatic = Chem.MolFromSmiles('c1ccccc1C')
        saturated = Chem.MolFromSmiles('C1CCCCC1C')
        self.assertTrue(fixed_substructure_match(aromatic, aromatic, list(range(6))))
        self.assertFalse(fixed_substructure_match(saturated, aromatic, list(range(6))))
        self.assertEqual(expected_fixed(aromatic), [6])

    def test_unique_float32_frame_conversion(self):
        rng = np.random.default_rng(23)
        x = rng.normal(size=(100, 3)).astype(np.float32)
        origin = np.array([-.06740643829107285, 28.96502113342285, 8.85383129119873], dtype=np.float32)
        world = x + origin
        recovered = recover_float32_origin(x, world, np.ones(len(x), bool))
        np.testing.assert_array_equal(recovered, origin)
        np.testing.assert_array_equal(x + recovered, world)

    def test_observation_is_not_causal_onset(self):
        result = suffix([{'value': False}, {'value': True}, {'value': False}, {'value': True}])
        self.assertIsNone(result['causal_onset'])
        self.assertTrue(result['boundary_is_confirmed_change'])

    def test_path_ids_cannot_escape_output(self):
        for value in ('../sample', '/tmp/sample', '.', '..', 'a/b'):
            with self.assertRaises(ValueError): safe_id(value)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory); (p/'keep').write_text('unchanged')
            with self.assertRaises(ValueError): check([], p)
            self.assertEqual((p/'keep').read_text(), 'unchanged')

    def test_missing_input_stops_before_experiment(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'manifest.json'
            p.write_text(json.dumps({'schema': 'scaffold_repair_manifest_v1', 'samples': [
                {'sdf': 'missing.sdf', 'reference_sdf': 'missing_ref.sdf', 'protein_pdb': 'missing.pdb'}]}))
            with self.assertRaises(ValueError): load_manifest(p)

    def test_no_flowr_runtime_imports(self):
        root = Path(__file__).resolve().parents[1] / 'scaffold_repair'
        for path in root.rglob('*.py'):
            for node in ast.walk(ast.parse(path.read_text())):
                modules = [a.name for a in node.names] if isinstance(node, ast.Import) else [node.module or ''] if isinstance(node, ast.ImportFrom) else []
                self.assertFalse(any(m == 'flowr' or m.startswith('flowr.') for m in modules), str(path))

    def test_geometry_source_remains_byte_identical(self):
        root = Path(__file__).resolve().parents[1]
        provenance = json.loads((root/'SOURCE_PROVENANCE.json').read_text())
        for item in provenance:
            if item.get('change') == 'byte-for-byte copy':
                self.assertEqual(hashlib.sha256((root/item['destination']).read_bytes()).hexdigest(), item['sha256'])


if __name__ == '__main__':
    unittest.main()

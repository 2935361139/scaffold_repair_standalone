"""Explicit evaluation scope; missing input is never a protein-free opt-in."""
from copy import deepcopy
from .io import require, resolve

ENVIRONMENT_FUNCTIONS = frozenset(('intermolecular_distance', 'volume_overlap'))
MOLECULAR_SCOPE = 'molecule_only'
PROTEIN_SCOPE = 'molecule_and_environment'


def protein_input(mode, value, base):
    require(mode in ('required', 'absent'), 'protein_mode must be required or absent')
    if mode == 'absent':
        require(not value, 'protein_mode=absent conflicts with a supplied PDB path')
        return None
    require(bool(value), 'Protein is required. Supply protein_pdb, or explicitly set protein_mode=absent')
    path = resolve(value, base)
    require(path.is_file(), f'Input file missing: {path}')
    return path


def molecular_config(config):
    """Subset of the original dock profile; every retained threshold is unchanged."""
    result = deepcopy(config)
    result['modules'] = [m for m in result['modules'] if m['function'] not in ENVIRONMENT_FUNCTIONS]
    for module in result['modules']:
        if module['function'] == 'loading':
            module['chosen_binary_test_output'] = [k for k in module['chosen_binary_test_output'] if k != 'mol_cond_loaded']
    return result


def omitted_checks(config):
    keys = []
    for module in config['modules']:
        for key in module.get('chosen_binary_test_output', []):
            if module['function'] in ENVIRONMENT_FUNCTIONS or key == 'mol_cond_loaded':
                keys.append(module.get('rename_outputs', {}).get(key, key).lower().replace(' ', '_'))
    return keys


def scope_record(protein, config):
    absent = protein is None
    skipped = {k: dict(status='not_evaluated', reason='no_protein_provided', passed=None)
               for k in omitted_checks(config)} if absent else {}
    return dict(evaluation_scope=MOLECULAR_SCOPE if absent else PROTEIN_SCOPE,
                protein_context_status='not_evaluated_no_protein' if absent else 'provided',
                not_evaluated_checks=skipped,
                success_meaning='Pass all evaluated molecular checks; protein environment NOT evaluated' if absent
                    else 'Pass the original molecular and protein-environment endpoint')


def scope_counts(rows, intervention=False):
    result = {}
    for row in rows:
        scope = row['evaluation_scope']
        bucket = result.setdefault(scope, dict(total=0, passed=0, failed=0))
        passed = row['selected_arm'] is not None if intervention else row['success']
        bucket['total'] += 1
        bucket['passed' if passed else 'failed'] += 1
        if intervention:
            bucket['rescued'] = bucket.get('rescued', 0) + int(passed and not row['original_success'])
    return result

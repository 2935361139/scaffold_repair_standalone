"""Snapshot inputs, check once, intervene once, retain audit evidence."""
import platform
import shutil
import traceback
import zipfile
from pathlib import Path
import numpy as np
import scipy
import rdkit
from rdkit import Chem
from posebusters import PoseBusters
from posebusters.tools.loading import safe_load_mol
from .io import require, read_json, write_json, sha, save_csv, safe_id
from .evaluation import evaluate, compare_endpoint, read_mol, save_mol, restore_flowr_decoder_representation
from .diagnosis import diagnose
from .geometry_intervention import SETTINGS, successful, route_failure, select_result
from .engine import ARMS, run_bank, plan_bank
from .receptor import prepare as prepare_receptor
from .topology import LIMITS
from .conformers import SETTINGS as CONFORMER_SETTINGS
from .coverage import assess, batch_metrics, CATALOG, SOURCES
from .context import scope_counts, molecular_config, MOLECULAR_SCOPE, PROTEIN_SCOPE


def versions():
    import pandas, yaml
    return dict(python=platform.python_version(), rdkit=rdkit.__version__, numpy=np.__version__,
                scipy=scipy.__version__, pandas=pandas.__version__, pyyaml=yaml.__version__)


def method_record():
    root = Path(__file__).resolve().parent.parent
    hashes = {str(p.relative_to(root)): sha(p) for name in ('scaffold_repair', 'posebusters')
              for p in sorted((root / name).rglob('*')) if p.suffix in ('.py', '.yml', '.yaml')}
    return dict(schema='standalone_scaffold_method_v2', versions=versions(), code_sha256=hashes,
                evaluator_config=PoseBusters('dock').config,
                molecular_evaluator_config=molecular_config(PoseBusters('dock').config),
                scope_policy='Explicit absent protein only; omitted environment checks are not evaluated, never passed',
                optimizer_settings=SETTINGS,
                topology_settings=LIMITS.copy(),
                conformer_settings=CONFORMER_SETTINGS.copy(),
                arms=ARMS, primary_priority=['A', 'B', 'AB', 'C'],
                new_generations=0, neural_interventions=0,
                stop_rule='One fixed bank per eligible sample; no retries or outcome-based parameter tuning',
                interpretation='V2: corrected receptor selections, bounded fixed-template/connectivity edits plus geometry controls. Same five endpoint conditions; new baselines required. No unique upstream cause claim.')


def snapshot(item, destination, receptor_options=None):
    destination.mkdir(parents=True)
    (destination / 'raw.sdf').write_bytes(item['raw_sdf'])
    save_mol(destination / 'reference.sdf', item['ref'])
    mode = item.get('protein_mode', 'required')
    require((mode == 'absent') == (item['protein'] is None), 'Protein mode/context mismatch')
    if item['protein'] is not None:
        shutil.copyfile(item['protein'], destination / 'protein.original.pdb')
        prepare_receptor(item['protein'], item['ref'], destination / 'protein.pdb', receptor_options)
    else:
        require(not receptor_options, 'RECEPTOR_OPTIONS cannot be applied without a protein')
        write_json(destination / 'protein_context.json', dict(status='not_evaluated', reason='explicit_no_protein'))
    np.savez_compressed(destination / 'coordinates.npz', xyz=item['xyz'], target=item['target'])
    if item['trajectory']:
        np.savez_compressed(destination / 'trajectory.npz', **item['trajectory'])
    meta = {k: item[k] for k in ('id', 'input_id', 'seed', 'slots', 'task', 'expected',
                                'trajectory_meta', 'imported_observations', 'provenance', 'raw_fixed_evidence')}
    meta['protein_mode'] = mode
    write_json(destination / 'sample.json', meta)


def restore(destination):
    item = read_json(destination / 'sample.json')
    with np.load(destination / 'coordinates.npz', allow_pickle=False) as arrays:
        item['xyz'], item['target'] = arrays['xyz'], arrays['target']
    item.update(mol=read_mol(destination / 'raw.sdf'), ref=read_mol(destination / 'reference.sdf'),
                protein=None if item.get('protein_mode') == 'absent' else destination / 'protein.pdb',
                raw_sdf=(destination / 'raw.sdf').read_bytes(), trajectory={})
    if item['provenance']['adapter'] == 'flowr_records':
        item['mol'] = restore_flowr_decoder_representation(item['mol'])
    Chem.SanitizeMol(item['ref'])
    if (destination / 'trajectory.npz').is_file():
        with np.load(destination / 'trajectory.npz', allow_pickle=False) as arrays:
            item['trajectory'] = {k: arrays[k] for k in arrays.files}
    return item


def environment(path, config, cache):
    if path is None:
        return None
    digest = sha(path)
    if digest not in cache:
        cache[digest] = safe_load_mol(path, **config.get('loading', {}).get('mol_cond', {}))
        require(cache[digest] is not None and cache[digest].GetNumConformers() == 1,
                f'Protein environment could not be loaded: {path}')
        require(np.isfinite(cache[digest].GetConformer().GetPositions()).all(), 'Nonfinite protein coordinates')
    return cache[digest]


def bundle(output):
    with zipfile.ZipFile(output / 'review.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(output.rglob('*')):
            if path.is_file() and path.name != 'review.zip' and not path.name.endswith('.tmp'):
                archive.write(path, path.relative_to(output))


def assert_sources_unchanged(items):
    for item in items:
        for path, digest in item['provenance']['source_hashes'].items():
            require(Path(path).is_file() and sha(path) == digest, f'Source changed during processing: {path}')


def check(items, output, receptor_options=None):
    output = Path(output).resolve()
    require(not output.exists(), f'Output already exists: {output}. Choose a new OUTPUT_DIR; existing results are preserved.')
    require(items and len({item['id'] for item in items}) == len(items), 'Empty or duplicated inputs')
    # A passing control must genuinely belong to the same reference/protein condition.
    conditions = {}
    for item in items:
        signature = (sha(item['protein']) if item['protein'] is not None else 'explicit_no_protein',
                     Chem.MolToMolBlock(item['ref']), tuple(item['task']['fixed_reference_indices']))
        previous = conditions.setdefault(item['input_id'], signature)
        require(previous == signature, f"input_id {item['input_id']} identifies inconsistent conditions")
    output.mkdir(parents=True)
    for name in ('success', 'failed', 'samples'):
        (output / name).mkdir()
    method = method_record()
    write_json(output / 'method.json', method)
    write_json(output / 'coverage_catalog.json', dict(categories=CATALOG,sources=SOURCES))
    rows, errors, cache = [], [], {}
    try:
        for number, item in enumerate(items, 1):
            rid = safe_id(item['id'])
            print(f'CHECK {number}/{len(items)} {rid}', flush=True)
            folder = output / 'samples' / rid
            snapshot(item, folder, receptor_options)
            current = restore(folder)
            row, mol = evaluate(current['mol'], current['xyz'][current['slots']], current['ref'],
                                current['task'], current['slots'], current['protein'], current['raw_fixed_evidence'])
            if item['expected'] is not None:
                # V1 environmental bugs are deliberately not reproduced. Preserve
                # chemistry checks and report legacy differences separately.
                for key in ('rdkit_valid','fc_valid','root_fixed_substructure_match','fixed_chemistry_preserved'):
                    require(row.get(key) == item['expected'].get(key), f'Imported chemical baseline mismatch: {key}')
                write_json(folder / 'legacy_evaluation_comparison.json', dict(previous=item['expected'],
                    changed_fields=[k for k in item['expected'] if row.get(k) != item['expected'][k]],
                    explanation='V2 fixes receptor selection; old and new environmental labels are not silently pooled.'))
            env = environment(current['protein'], method['evaluator_config'], cache)
            diagnostic = diagnose(current, mol, row, env, method['evaluator_config'])
            write_json(folder / 'evaluation.json', row)
            write_json(folder / 'diagnosis.json', diagnostic)
            write_json(folder / 'research_assessment.json',assess(current,mol,row,method['evaluator_config']))
            route = route_failure(row)
            result = dict(id=rid, input_id=item['input_id'], seed=item['seed'], success=successful(row),
                          evaluation_scope=row['evaluation_scope'], success_meaning=row['success_meaning'],
                          route=route, evaluation=row, localized_pair_count=len(diagnostic['pairs']),
                          timing_status=diagnostic['temporal']['status'])
            rows.append(result)
            classified = output / ('success' if result['success'] else 'failed') / rid
            classified.mkdir()
            save_mol(classified / 'molecule.sdf', mol)
            if mol is None:
                (classified / 'molecule.sdf').write_bytes(item['raw_sdf'])
            write_json(classified / 'sample.json', result)
        assert_sources_unchanged(items)
    except Exception as exc:
        errors.append(dict(error=str(exc), traceback=traceback.format_exc(), sample=items[len(rows)]['id'] if len(rows) < len(items) else None))
    summary = dict(completed=not errors and len(rows) == len(items), requested_samples=len(items),
        evaluated_samples=len(rows), passed=sum(r['success'] for r in rows),
        failed=sum(not r['success'] for r in rows),
        supported_geometry_failures=sum(r['route']['status'] == 'supported_geometry' for r in rows),
        supported_topology_failures=sum(r['route']['status'] == 'supported_topology' for r in rows),
        unsupported_failures=sum(r['route']['status'] == 'unsupported' for r in rows),
        errors=errors, samples=rows, new_generations=0,
        by_evaluation_scope=scope_counts(rows),
        interpretation='Pass/fail applies only to evaluation_scope. Molecule-only passes do not establish protein compatibility; compare scopes separately.')
    write_json(output / 'summary.json', summary)
    write_json(output / 'batch_metrics.json',batch_metrics([r['evaluation']['molecule_smiles'] for r in rows],len(items)))
    save_csv(output / 'samples.csv', [dict(id=r['id'], input_id=r['input_id'], success=r['success'],
        evaluation_scope=r['evaluation_scope'], not_evaluated_checks=list(r['evaluation']['not_evaluated_checks']),
        route=r['route']['status'], failed_checks=[k for k in ('rdkit_valid','fc_valid','root_fixed_substructure_match','fixed_chemistry_preserved') if r['evaluation'].get(k) is not True] + r['evaluation']['posebusters_failed_checks'],
        localized_pairs=r['localized_pair_count'], timing=r['timing_status'],
        murcko_scaffold_changed=r['evaluation']['murcko_scaffold_changed']) for r in rows])
    if summary['completed']:
        hashes = {str(path.relative_to(output)): sha(path) for path in sorted(output.rglob('*')) if path.is_file()}
        write_json(output / 'integrity.json', hashes)
    bundle(output)
    require(summary['completed'], f"Evaluation stopped: {errors[0]['error'] if errors else 'incomplete'}; see {output / 'review.zip'}")
    return summary


def verify_checked(root):
    require(read_json(root / 'summary.json')['completed'] is True, 'Incomplete checked data cannot enter interventions')
    hashes = read_json(root / 'integrity.json')
    for relative, digest in hashes.items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root) and path.is_file() and sha(path) == digest,
                f'Checked evidence changed: {relative}. Run a fresh check in a new directory.')
    expected = read_json(root / 'method.json')
    current = method_record()
    for key in ('schema', 'code_sha256', 'evaluator_config', 'optimizer_settings', 'topology_settings', 'conformer_settings', 'arms', 'versions'):
        require(current[key] == expected[key], f'Method/environment changed between stages: {key}')


def intervene(failed_dir, output, failure_ids=None):
    failed_dir, output = Path(failed_dir).resolve(), Path(output).resolve()
    root = failed_dir.parent
    require(failed_dir.name == 'failed' and failed_dir.is_dir() and (root / 'integrity.json').is_file(),
            'INTERVENTION_INPUT_DIR must be the failed directory produced by this standalone tool')
    require(not output.exists(), f'Intervention output exists: {output}. Preserve it and choose a new OUTPUT_DIR.')
    verify_checked(root)
    checked = read_json(root / 'summary.json')['samples']
    available = {r['id'] for r in checked if not r['success']}
    requested = available if failure_ids is None else set(failure_ids)
    require(failure_ids is None or len(failure_ids) == len(requested), 'Duplicate requested failure IDs')
    require(requested <= available, 'A requested failure is absent or originally passed')
    active_inputs = {r['input_id'] for r in checked if r['id'] in requested and r['route']['status'] in ('supported_geometry','supported_topology')}
    active_ids = {r['id'] for r in checked if r['input_id'] in active_inputs and
                  (r['id'] in requested or r['success']) and r['route']['status'] != 'unsupported'}
    output.mkdir(parents=True)
    for name in ('success', 'failed', 'arms'):
        (output / name).mkdir()
    method = method_record()
    write_json(output / 'method.json', method)
    write_json(output / 'coverage_catalog.json', dict(categories=CATALOG,sources=SOURCES))
    records, errors, conclusions, cache, fatal = [], [], [], {}, []
    plans, planned = {}, 0
    try:
        # Reproduce all checked baselines before any treated branch.
        items = {}
        for r in checked:
            item = restore(root / 'samples' / r['id'])
            baseline, _ = evaluate(item['mol'], item['xyz'][item['slots']], item['ref'], item['task'], item['slots'], item['protein'], item['raw_fixed_evidence'])
            require(baseline == r['evaluation'], f"Baseline no longer reproduces: {r['id']}")
            items[r['id']] = item
        for r in checked:
            if r['id'] in active_ids:
                plans[r['id']] = plan_bank(items[r['id']], r['evaluation'])
        planned = sum(len(p['arm_names']) for p in plans.values())
        write_json(output / 'selection.json', dict(source_generated=str(root), source_integrity_sha256=sha(root / 'integrity.json'),
            requested_failures=sorted(requested), active_runs=sorted(active_ids), planned_arms=planned,
            per_sample={rid:dict(arms=p['arm_names'],topology=p['audit']) for rid,p in plans.items()}))
        for r in checked:
            rid = r['id']
            item = items[rid]
            if rid in active_ids:
                env = environment(item['protein'], method['evaluator_config'], cache)
                branches, branch_errors, conclusion = run_bank(item, r['evaluation'], env, method['evaluator_config'], output / 'arms' / rid, plans[rid])
                records.extend(branches)
                errors.extend(branch_errors)
            else:
                conclusion = dict(run_id=rid, **select_result(r['evaluation'], {}))
            conclusion.update(original_success=r['success'], requested_failure=rid in requested,
                              evaluation_scope=r['evaluation_scope'], success_meaning=r['success_meaning'],
                              not_evaluated_checks=r['evaluation']['not_evaluated_checks'],
                              intervention_executed=rid in plans and len(plans[rid]['arm_names']) > 0, route=r['route']['status'])
            if not r['success'] and rid not in requested:
                conclusion['result'] = 'not_requested'
            elif r['route']['status'] == 'unsupported':
                conclusion['result'] = 'unsupported_failure'
            conclusions.append(conclusion)
            selected = conclusion['selected_arm']
            dest = output / ('success' if selected else 'failed') / rid
            dest.mkdir()
            if selected and selected != 'baseline':
                shutil.copyfile(output / 'arms' / rid / (selected + '.sdf'), dest / 'molecule.sdf')
                np.save(dest / 'coordinates.npy', np.load(output / 'arms' / rid / (selected + '.npy'), allow_pickle=False)[item['slots']], allow_pickle=False)
            else:
                shutil.copyfile(root / ('success' if r['success'] else 'failed') / rid / 'molecule.sdf', dest / 'molecule.sdf')
                np.save(dest / 'coordinates.npy', item['xyz'][item['slots']], allow_pickle=False)
            write_json(dest / 'sample.json', conclusion)
            # Assess the selected full-precision structure, never just rounded SDF coordinates.
            from .evaluation import with_coords
            chosen_mol=read_mol(dest/'molecule.sdf')
            if chosen_mol is not None:
                chosen_mol=with_coords(chosen_mol,np.load(dest/'coordinates.npy',allow_pickle=False))
            chosen_row=next((v['evaluation'] for v in records if v['run_id']==rid and v['arm']==selected),r['evaluation'])
            write_json(dest/'research_assessment.json',assess(item,chosen_mol,chosen_row,method['evaluator_config']))
        verify_checked(root)
    except Exception as exc:
        fatal.append(dict(error=str(exc), traceback=traceback.format_exc()))
    rescued = [r for r in conclusions if not r['original_success'] and r['selected_arm'] is not None]
    summary = dict(completed=not fatal and len(conclusions) == len(checked),
        total_samples=len(checked), original_passed=sum(r['success'] for r in checked),
        original_failed=len(available), requested_failures=len(requested),
        rescued=len(rescued), rescued_without_relaxing_fixed_coordinates=sum(not r['fixed_coordinate_relaxation'] for r in rescued),
        rescued_with_relaxed_fixed_coordinates=sum(r['fixed_coordinate_relaxation'] for r in rescued),
        rescued_with_graph_edits=sum(bool(r.get('graph_changed')) for r in rescued),
        attempted_failures=sum(r['requested_failure'] and r['intervention_executed'] for r in conclusions),
        final_passed=sum(r['selected_arm'] is not None for r in conclusions),
        final_failed=sum(r['selected_arm'] is None for r in conclusions),
        unsupported_failures=sum(r['route']['status'] == 'unsupported' for r in checked),
        planned_arms=planned, completed_arms=len(records), branch_errors=len(errors),
        all_planned_arms_evaluated=len(records) == planned and not errors,
        fatal_errors=fatal, new_generations=0, neural_interventions=0,
        by_evaluation_scope=scope_counts(conclusions, intervention=True),
        interpretation='V2 whole-endpoint repairs. Graph edits and relaxed fixed coordinates are reported separately; unsupported, unattempted, evaluated failures and runtime errors are distinct.')
    write_json(output / 'summary.json', summary)
    write_json(output / 'results.json', dict(summary=summary, records=records, conclusions=conclusions, branch_errors=errors))
    chosen_smiles=[]
    for c in conclusions:
        original=next(r for r in checked if r['id']==c['run_id'])
        chosen=next((r['evaluation'] for r in records if r['run_id']==c['run_id'] and r['arm']==c['selected_arm']),original['evaluation'])
        chosen_smiles.append(chosen['molecule_smiles'])
    write_json(output/'batch_metrics.json',batch_metrics(chosen_smiles,len(checked)))
    save_csv(output / 'samples.csv', conclusions)
    save_csv(output / 'arms.csv', [dict(id=r['run_id'], arm=r['arm'], passed=r['full_success'],
        evaluation_scope=r['evaluation']['evaluation_scope'],
        improved=r['checks_improved'], worsened=r['checks_worsened'], failed_checks=r['evaluation']['posebusters_failed_checks']) for r in records])
    bundle(output)
    require(summary['completed'], f"Intervention stopped: {fatal[0]['error'] if fatal else 'incomplete'}; see {output / 'review.zip'}")
    return summary

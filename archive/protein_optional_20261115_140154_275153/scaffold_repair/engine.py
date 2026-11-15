"""Predeclared graph, anchored-conformer and original geometry intervention bank."""
import hashlib
import traceback
import numpy as np
from .io import require, write_json
from .evaluation import evaluate, with_coords, save_mol
from .legacy_functions import coordinate_arm
from .geometry_intervention import environment_terms, relax_generated, successful, select_result, classify_rescue
from .diagnosis import clash_pairs

ARMS = ['A', 'B', 'AB', 'C', 'UFF_only', 'A_UFF_only', 'C_UFF_only',
        'A_random_0', 'A_random_1', 'A_random_2', 'B_random_0', 'B_random_1', 'B_random_2',
        'AB_random_0', 'AB_random_1', 'AB_random_2', 'C_random_0', 'C_random_1', 'C_random_2']


def run_geometry_bank(item, baseline, environment, config, destination):
    destination.mkdir(parents=True, exist_ok=True)
    rid = item['id']
    xyz0, target, slots = item['xyz'], item['target'], item['slots']
    fixed = item['task']['fixed_sampler_indices']
    generated = [s for s in slots if s not in fixed]
    fixed_heavy = [i for i, s in enumerate(slots) if s in fixed]
    require(generated and len(fixed_heavy) == len(fixed), 'Invalid fixed/generated mapping')
    baseline_row, base_mol = evaluate(item['mol'], xyz0[slots], item['ref'], item['task'], slots, item['protein'], item['raw_fixed_evidence'])
    require(baseline_row == baseline, 'Baseline changed between checking and intervention')
    terms = environment_terms(base_mol, environment, config)
    xyzA, receiptA = coordinate_arm(xyz0, target, fixed, 'restore_fixed', 0)
    designs = {'A': (xyzA, receiptA)}
    errors, records, evaluated = [], [], {}

    def error(arm, exc):
        result = dict(run_id=rid, arm=arm, status='not_evaluated', error=str(exc), traceback=traceback.format_exc())
        errors.append(result)
        write_json(destination / (arm + '_error.json'), result)

    for name, start, use_env, move_fixed in [('B', xyz0, True, False), ('AB', xyzA, True, False),
            ('UFF_only', xyz0, False, False), ('A_UFF_only', xyzA, False, False),
            ('C', xyz0, True, True), ('C_UFF_only', xyz0, False, True)]:
        print(f'OPTIMIZE {rid} {name}', flush=True)
        try:
            mol = with_coords(base_mol, start[slots])
            heavy_xyz, receipt = relax_generated(mol, [] if move_fixed else fixed_heavy, terms, use_environment=use_env)
            result = start.copy()
            for i, s in enumerate(slots):
                if move_fixed or s in generated:
                    result[s] = heavy_xyz[i]
            receipt.update(base='A' if name in ('AB', 'A_UFF_only') else 'baseline',
                fixed_coordinate_constraint_relaxed=move_fixed, condition_fixed_heavy_indices=fixed_heavy,
                moved_heavy_slots=slots if move_fixed else generated,
                max_fixed_change_from_baseline_A=float(np.linalg.norm(result[fixed].astype(float) - xyz0[fixed].astype(float), axis=1).max()))
            if move_fixed:
                receipt['only_generated_heavy_coordinates_changed'] = False
                receipt['optimized_heavy_indices'] = receipt.pop('generated_heavy_indices')
                receipt['coordinate_scope'] = 'All ligand heavy atoms, including condition-fixed atoms'
            designs[name] = result, receipt
        except Exception as exc:
            error(name, exc)
    for name in ('A', 'B', 'AB', 'C'):
        if name not in designs:
            continue
        start = xyzA if name == 'AB' else xyz0
        selection = fixed if name == 'A' else (slots if name == 'C' else generated)
        for number in range(3):
            label = f'{name}_random_{number}'
            seed = int.from_bytes(hashlib.sha256(('20261029:' + rid + ':' + label).encode()).digest()[:8], 'little')
            try:
                result, receipt = coordinate_arm(start, designs[name][0], selection, 'random_direction', seed)
                receipt.update(matched_arm=name, base='A' if name == 'AB' else 'baseline')
                designs[label] = result, receipt
            except Exception as exc:
                error(label, exc)
    for arm in ARMS:
        if arm not in designs:
            if not any(e['arm'] == arm for e in errors):
                error(arm, 'Parent intervention unavailable')
            continue
        try:
            xyz, receipt = designs[arm]
            row, mol = evaluate(base_mol, xyz[slots], item['ref'], item['task'], slots, item['protein'])
            for key in ('rdkit_valid', 'fc_valid', 'root_fixed_substructure_match', 'fixed_chemistry_preserved',
                        'murcko_scaffold_changed', 'molecule_smiles', 'scaffold_smiles'):
                require(row[key] == baseline[key], f'Coordinate-only arm changed chemistry: {key}')
            pairs = clash_pairs(mol, item['ref'], item['task'], slots, environment, config, row)
            old, new = baseline['posebusters_checks'], row['posebusters_checks']
            record = dict(run_id=rid, input_id=item['input_id'], seed=item['seed'], arm=arm, receipt=receipt,
                evaluation=row, pairs=pairs, full_success=successful(row),
                operational_rescue=not successful(baseline) and successful(row),
                checks_improved=[k for k in old if old[k] is False and new[k] is True],
                checks_worsened=[k for k in old if old[k] is True and new[k] is False])
            save_mol(destination / (arm + '.sdf'), mol)
            np.save(destination / (arm + '.npy'), xyz, allow_pickle=False)
            write_json(destination / (arm + '.json'), record)
            records.append(record)
            evaluated[arm] = row
            print(f'ARM {rid} {arm}: success={successful(row)}', flush=True)
        except Exception as exc:
            error(arm, exc)
    conclusion = dict(run_id=rid, **classify_rescue(baseline, evaluated), **select_result(baseline, evaluated))
    write_json(destination / 'conclusion.json', conclusion)
    return records, errors, conclusion


def plan_bank(item, baseline):
    from .topology import plan
    from .conformers import propose
    needs_graph = not all(baseline.get(k) is True for k in
        ('rdkit_valid','fc_valid','root_fixed_substructure_match','fixed_chemistry_preserved'))
    candidates, audit = plan(item) if needs_graph else ([], dict(status='graph_already_passes', candidates=[]))
    rebuilt, embedding_audit = propose(item,candidates) if not successful(baseline) else ([],dict(status='baseline_already_passes'))
    candidates += rebuilt
    audit['conformer_proposals']=embedding_audit
    geometry = baseline['rdkit_valid'] is True
    names = list(ARMS) if geometry else []
    for candidate in candidates:
        names += [candidate['name']+suffix for suffix in ('','_B','_C','_UFF_control')]
    return dict(candidates=candidates, audit=audit, geometry=geometry, arm_names=names)


def run_bank(item, baseline, environment, config, destination, plan=None):
    from .topology import graph_edits, graph_key
    plan = plan or plan_bank(item, baseline)
    destination.mkdir(parents=True, exist_ok=True)
    write_json(destination / 'planning.json', dict(**plan['audit'], planned_arms=plan['arm_names']))
    records, errors = [], []
    if plan['geometry']:
        # Original geometry bank serves both as a repair and as the graph-unchanged control.
        records, errors, _ = run_geometry_bank(item, baseline, environment, config, destination)
        for record in records:
            record.update(is_primary=record['arm'] in ('A','B','AB','C'), graph_intervention=None,
                          fixed_coordinate_relaxation=record['arm'] in ('C','C_UFF_only') or record['arm'].startswith('C_random_'),
                          graph_edits=dict(bonds=[],charges=[]))
    xyz0, slots = item['xyz'], item['slots']
    fixed_heavy = [i for i,s in enumerate(slots) if s in item['task']['fixed_sampler_indices']]
    for proposal in plan['candidates']:
        start = proposal.get('coordinates',xyz0[slots])
        base = with_coords(proposal['mol'], start)
        for suffix in ('', '_B', '_C', '_UFF_control'):
            arm = proposal['name'] + suffix
            try:
                coords = start.copy()
                receipt = dict(graph_edits=proposal['edits'], family=proposal['family'],
                               geometry='none', fixed_coordinate_constraint_relaxed=suffix == '_C')
                if suffix:
                    terms = environment_terms(base, environment, config)
                    coords, optim = relax_generated(base, [] if suffix == '_C' else fixed_heavy,
                        terms, use_environment=suffix != '_UFF_control')
                    receipt.update(optimization=optim, geometry=suffix[1:])
                xyz = xyz0.copy()
                xyz[slots] = coords
                row, mol = evaluate(base, coords, item['ref'], item['task'], slots, item['protein'])
                require(graph_key(mol) == graph_key(base), 'Coordinate optimization changed the proposed graph')
                before = {k:v for k,v in baseline.items() if isinstance(v,bool)}
                before.update(baseline.get('posebusters_checks') or {})
                after = {k:v for k,v in row.items() if isinstance(v,bool)}
                after.update(row.get('posebusters_checks') or {})
                changed_graph=bool(proposal['edits']['bonds'] or proposal['edits']['charges'])
                record = dict(run_id=item['id'], input_id=item['input_id'], seed=item['seed'], arm=arm,
                    is_primary=suffix != '_UFF_control', graph_intervention=proposal['family'] if changed_graph else None,
                    conformer_reconstruction=proposal['family']=='G',proposal_receipt=proposal.get('receipt'),
                    graph_edits=graph_edits(item['mol'],mol), fixed_coordinate_relaxation=suffix=='_C',
                    receipt=receipt, evaluation=row, full_success=successful(row),
                    operational_rescue=not successful(baseline) and successful(row),
                    checks_improved=[k for k,v in before.items() if v is False and after.get(k) is True],
                    checks_worsened=[k for k,v in before.items() if v is True and after.get(k) is False],
                    pairs=clash_pairs(mol,item['ref'],item['task'],slots,environment,config,row))
                save_mol(destination / (arm+'.sdf'),mol)
                np.save(destination / (arm+'.npy'),xyz,allow_pickle=False)
                write_json(destination / (arm+'.json'),record)
                records.append(record)
                print(f'ARM {item["id"]} {arm}: success={successful(row)}',flush=True)
            except Exception as exc:
                entry=dict(run_id=item['id'],arm=arm,status='not_evaluated',error=str(exc),traceback=traceback.format_exc())
                errors.append(entry)
                write_json(destination / (arm+'_error.json'),entry)
    # Persist normalized roles for old geometry branches as well.
    for record in records:
        write_json(destination / (record['arm']+'.json'),record)
    primary = [r for r in records if r['is_primary'] and r['full_success']]
    primary.sort(key=lambda r:(r['fixed_coordinate_relaxation'],len(r['graph_edits']['bonds'])+len(r['graph_edits']['charges']),
                               plan['arm_names'].index(r['arm'])))
    selected = 'baseline' if successful(baseline) else (primary[0]['arm'] if primary else None)
    conclusion = dict(run_id=item['id'],selected_arm=selected,
        result='original_pass' if selected=='baseline' else ('rescued_v2' if selected else 'unresolved_within_budget'),
        all_successful_interventions=[r['arm'] for r in primary],
        fixed_coordinate_relaxation=bool(selected and selected!='baseline' and primary[0]['fixed_coordinate_relaxation']),
        graph_changed=bool(selected and selected!='baseline' and primary[0]['graph_intervention']),
        successful_mechanisms=[dict(arm=r['arm'],graph=r['graph_intervention'],
            fixed_coordinates_relaxed=r['fixed_coordinate_relaxation']) for r in primary],
        planned_arms=len(plan['arm_names']),evaluated_arms=len(records),branch_errors=len(errors),
        candidate_status=plan['audit']['status'],
        interpretation='Only whole-endpoint passes count. Combined graph/geometry rescue does not establish either factor alone or a unique neural cause.')
    write_json(destination / 'conclusion.json',conclusion)
    return records,errors,conclusion

"""Finite endpoint localization and optional observed trajectory timing."""
import numpy as np
from rdkit import Chem
from posebusters.modules.intermolecular_distance import check_intermolecular_distance
from .io import require
from .regions import atom_regions, fixed_pair_comparison
from .legacy_functions import series, suffix


def environment_atom(atom):
    info = atom.GetPDBResidueInfo()
    return dict(environment_atom_index=atom.GetIdx(), environment_element=atom.GetSymbol(),
        pdb_serial=info.GetSerialNumber() if info else None,
        atom_name=info.GetName().strip() if info else None,
        residue=info.GetResidueName().strip() if info else None,
        residue_number=info.GetResidueNumber() if info else None,
        chain=info.GetChainId().strip() if info else None)


def clash_pairs(mol, ref, task, slots, environment, config, row):
    if mol is None or row['rdkit_valid'] is not True or row['posebusters_checks'] is None:
        return []
    mapping = dict(zip(task['fixed_sampler_indices'], task['fixed_reference_indices']))
    if not set(mapping) <= set(slots):
        return []
    anchors = {a['fixed_sampler_index'] for a in task['reference_attachment_sites']}
    edges = [(slots[b.GetBeginAtomIdx()], slots[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    roles = atom_regions(slots, edges, set(mapping), anchors)
    mol_h = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    ref_h = Chem.AddHs(Chem.Mol(ref), addCoords=True)
    pairs = []
    for module in config['modules']:
        if module['function'] != 'intermolecular_distance':
            continue
        label = module['rename_outputs'].get('no_clashes', 'no_clashes').lower().replace(' ', '_')
        kwargs = module.get('parameters', {})
        result = check_intermolecular_distance(mol_h, environment, **kwargs)
        require(bool(result['results']['no_clashes']) == row['posebusters_checks'][label], 'Clash details disagree with evaluator')
        reference = check_intermolecular_distance(ref_h, environment, **kwargs)
        ref_pairs = {(int(x['ligand_atom_id']), int(x['protein_atom_id']))
                     for x in reference['details'].to_dict('records') if bool(x['clash'])}
        for record in result['details'].to_dict('records'):
            if not bool(record['clash']):
                continue
            i, other = int(record['ligand_atom_id']), int(record['protein_atom_id'])
            require(i < len(slots), 'Clash unexpectedly assigned to auxiliary hydrogen')
            slot = slots[i]
            pairs.append(dict(check=label, ligand_sdf_atom_index=i, sampler_atom_index=slot,
                fixed_reference_atom_index=mapping.get(slot), element=mol.GetAtomWithIdx(i).GetSymbol(),
                **roles[slot], **environment_atom(environment.GetAtomWithIdx(other)),
                distance_A=float(record['distance']), relative_distance=float(record['relative_distance']),
                native_clash_cutoff=kwargs['clash_cutoff'],
                reference_pair_comparison=fixed_pair_comparison(slot, other, mapping, ref_pairs),
                causal_status='candidate_not_unique_upstream_cause'))
    return pairs


def temporal_evidence(item, pairs, environment, config):
    if not item['trajectory']:
        return dict(status='unavailable_final_SDF_only', observations=[],
                    explanation='No saved trajectory: generation step and neural module cannot be inferred from final SDF.')
    modules = {m['rename_outputs']['no_clashes'].lower().replace(' ', '_'): m.get('parameters', {})
               for m in config['modules'] if m['function'] == 'intermolecular_distance'}
    table = Chem.GetPeriodicTable()
    output = []
    for pair in pairs:
        params = modules[pair['check']]
        radius = table.GetRvdw if params.get('radius_type', 'vdw') == 'vdw' else table.GetRcovalent
        denom = (radius(pair['element']) + radius(pair['environment_element'])) * params.get('radius_scale', 1.)
        env_xyz = environment.GetConformer().GetPositions()[pair['environment_atom_index']]
        slot = pair['sampler_atom_index']
        for stream, metadata in item['trajectory_meta']['streams'].items():
            a = item['trajectory']
            points = series(a[stream + '_coords'], a[stream + '_atom_classes'], a[stream + '_valid'],
                            item['trajectory_meta']['atom_vocabulary'], slot, pair['element'],
                            env_xyz, denom, params['clash_cutoff'], metadata)
            output.append(dict(check=pair['check'], stream=stream, sampler_atom_index=slot,
                environment_atom_index=pair['environment_atom_index'], points=points, **suffix(points),
                internal_module=None, causal_status='not_tested',
                scope='Observed timing of this final clash pair; not a proven causal onset or earliest error anywhere.'))
    return dict(status='observed_saved_states', observations=output)


def diagnose(item, mol, row, environment, config):
    pairs = clash_pairs(mol, item['ref'], item['task'], item['slots'], environment, config, row)
    localized = {p['check'] for p in pairs}
    unlocalized = [dict(check=k, status='endpoint_failure_not_assigned_to_one_atom')
                   for k in row.get('posebusters_failed_checks', []) if k not in localized]
    for key in ('rdkit_valid', 'fc_valid', 'root_fixed_substructure_match', 'fixed_chemistry_preserved'):
        if row.get(key) is not True:
            unlocalized.append(dict(check=key, status='see_graph_diagnosis_and_bounded_topology_plan'))
    from .topology import bond_map, local_fixed, bridge_pairs
    graph = dict(status='unreadable' if mol is None else 'observed_final_graph',
                 components=[], fixed_bond_changes=[], proposed_cross_component_pairs=[])
    candidates = []
    if mol is not None:
        mapping = local_fixed(item)
        graph['components'] = [list(c) for c in Chem.GetMolFrags(mol)]
        generated_bonds, reference_bonds = bond_map(mol), bond_map(item['ref'])
        for i, j in row.get('fixed_internal_bond_changed_pairs', []):
            a, b = item['slots'].index(i), item['slots'].index(j)
            graph['fixed_bond_changes'].append(dict(generated_atoms=[a,b],
                reference_atoms=[mapping[a],mapping[b]],
                expected=reference_bonds.get(tuple(sorted((mapping[a],mapping[b])))),
                observed=generated_bonds.get(tuple(sorted((a,b))))))
        if graph['fixed_bond_changes'] or row.get('fixed_charge_changed_sites'):
            candidates.append(dict(kind='F', hypothesis='Mapped fixed-template bond/charge inconsistency',
                                   evidence=graph['fixed_bond_changes']))
        if len(graph['components']) > 1:
            graph['proposed_cross_component_pairs'] = [dict(atoms=[a,b], distance_A=d, distance_score=score)
                for score,a,b,d in bridge_pairs(mol,set(mapping))]
            candidates.append(dict(kind='L', hypothesis='Disconnected components may admit a local covalent connection',
                                   component_count=len(graph['components'])))
        if row.get('fixed_max_displacement_A',0) > 1e-5:
            candidates.append(dict(kind='A', hypothesis='Observed fixed-coordinate drift',
                                   max_displacement_A=row['fixed_max_displacement_A']))
    if any(k != 'all_atoms_connected' for k in row.get('posebusters_failed_checks', [])):
        candidates.append(dict(kind='geometry', hypothesis='Reported internal or environmental geometry failures',
                               failed_checks=row['posebusters_failed_checks']))
    return dict(pairs=pairs, unlocalized=unlocalized,
                graph=graph,
                temporal=temporal_evidence(item, pairs, environment, config),
                imported_observations=item.get('imported_observations'),
                candidates=candidates,
                interpretation='Observed graph/geometry failures select bounded intervention families. A successful repair is not proof of a unique upstream cause.')

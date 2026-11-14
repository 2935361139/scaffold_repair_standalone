"""Import completed records, never import FLOWR or load a model checkpoint."""
from pathlib import Path
import numpy as np
from rdkit import Chem
from ..io import read_json, require, sha, safe_id
from ..evaluation import read_mol, expected_fixed, restore_flowr_decoder_representation


def recover_float32_origin(centered, world, valid):
    """Recover a UNIQUE recorded float32 translation from exact rounding intervals.

    Refuse ambiguous translations. This is a representation conversion, not a fit,
    alignment or intervention. Every final coordinate must reproduce bit-for-bit.
    """
    require(centered.dtype == world.dtype == np.dtype('float32'), 'Origin recovery requires float32 records')
    x, y = centered[valid], world[valid]
    require(len(x) > 0 and np.isfinite(x).all() and np.isfinite(y).all(), 'Invalid terminal frame')
    origin = []
    for k in range(3):
        lo = float(np.max((y[:, k].astype(float) + np.nextafter(y[:, k], -np.inf).astype(float)) / 2 - x[:, k].astype(float)))
        hi = float(np.min((y[:, k].astype(float) + np.nextafter(y[:, k], np.inf).astype(float)) / 2 - x[:, k].astype(float)))
        require(lo <= hi, 'No consistent recorded origin; do not align approximately')
        candidate = np.float32(lo)
        if float(candidate) < lo:
            candidate = np.nextafter(candidate, np.float32(np.inf))
        matches = []
        for _ in range(4096):
            if float(candidate) > hi:
                break
            if np.array_equal(x[:, k] + candidate, y[:, k]):
                matches.append(candidate)
            if len(matches) > 1:
                break
            candidate = np.nextafter(candidate, np.float32(np.inf))
        else:
            raise ValueError('Origin recovery interval is too wide; provide original inputs.pt')
        require(len(matches) == 1, 'Recorded origin is ambiguous; retain original inputs.pt or use an explicit manifest')
        origin.append(matches[0])
    origin = np.asarray(origin, dtype=np.float32)
    require(np.array_equal(x + origin, y), 'Recovered origin does not reproduce final coordinates')
    return origin


def origin_for(folder, arrays):
    # Torch is optional and used only for reading an existing tensor file.
    inputs = folder / 'inputs.pt'
    if inputs.is_file() and not inputs.read_bytes()[:80].startswith(b'version https://git-lfs'):
        try:
            import torch
        except ImportError:
            pass
        else:
            record = torch.load(inputs, map_location='cpu', weights_only=True)
            require(len(record['pocket_com']) == 1, 'Expected one recorded origin')
            origin = record['pocket_com'][0].detach().cpu().numpy().reshape(3)
            valid = arrays['raw_final_valid']
            require(np.array_equal(arrays['sampler_state_coords'][-1][valid] + origin,
                                   arrays['raw_final_coords'][valid]), 'Recorded origin does not reproduce endpoint')
            return origin, 'original_inputs.pt', [inputs]
    return recover_float32_origin(arrays['sampler_state_coords'][-1], arrays['raw_final_coords'],
                                  arrays['raw_final_valid']), 'unique_float32_origin_verified_bitwise', []


def load(root, protein_override=None):
    root = Path(root).resolve()
    require(root.is_dir(), f'FLOWR generated directory not found: {root}')
    groups = sorted(p for p in root.glob('group_*') if p.is_dir()) if not (root / 'quality/quality_report.json').is_file() else [root]
    require(groups, 'No completed group directories found')
    items = []
    for group in groups:
        quality_path = group / 'quality/quality_report.json'
        require(quality_path.is_file(), f'Missing quality report: {quality_path}')
        quality = read_json(quality_path)
        require(quality.get('completed') is True and not quality.get('errors'), 'Original quality evaluation incomplete')
        diagnostics = {}
        study_path = group / 'localization_v2/study/diagnosis.json'
        study = {r['run_id']: r for r in read_json(study_path)['runs']} if study_path.is_file() else {}
        for path in sorted((group / 'diagnostics').glob('*/connection_report.json')):
            report = read_json(path)
            rid = Path(report['run']).name
            require(rid not in diagnostics, 'Duplicate diagnostic run ID')
            diagnostics[rid] = (path, report)
        rows = {Path(r['run_directory']).name: r for r in quality['rows'] if r['stage'] == 'raw_without_repair'}
        require(len(rows) == len(quality['runs']), 'Original sample count differs from endpoint rows')
        for info in quality['runs']:
            rid = safe_id(Path(info['run_directory']).name)
            folder = group / 'runs' / rid / 'research_record/attempt_000'
            require(rid in diagnostics and rid in rows, f'Missing evidence for {rid}')
            report_path, report = diagnostics[rid]
            array_path = report_path.with_name('connection_arrays.npz')
            with np.load(array_path, allow_pickle=False) as source:
                arrays = {k: source[k] for k in source.files}
            task_path = folder / 'task_definition.json'
            reference_path = folder / 'reference_preprocessed.sdf'
            sdf = folder / 'decoded_without_repair.sdf'
            for path, digest in ((task_path, info['task_sha256']), (reference_path, info['reference_sdf_sha256']),
                                 (sdf, rows[rid]['input_sdf_sha256'])):
                require(path.is_file() and sha(path) == digest, f'Original evidence missing or changed: {path}')
            task = read_json(task_path)
            require(task == report['task_definition'], 'Diagnostic task definition differs')
            ref = read_mol(reference_path)
            require(ref is not None, 'Unreadable reference')
            Chem.SanitizeMol(ref)
            require(sorted(task['fixed_reference_indices']) == expected_fixed(ref),
                    'FLOWR fixed mask differs from original scaffold-hopping definition')
            mol = restore_flowr_decoder_representation(read_mol(sdf))
            valid, classes = arrays['raw_final_valid'], arrays['raw_final_atom_classes']
            vocab = report['atom_vocabulary']
            slots = [i for i in range(len(valid)) if valid[i] and vocab[int(classes[i])] != '<NOATOM>']
            xyz = arrays['raw_final_coords']
            if mol is not None:
                require(mol.GetNumAtoms() == len(slots), 'SDF atom count differs from saved coordinates')
                require(all(a.GetSymbol() == vocab[int(classes[s])] for a, s in zip(mol.GetAtoms(), slots)), 'SDF element order mismatch')
                require(np.max(np.abs(mol.GetConformer().GetPositions() - xyz[slots])) <= .002, 'SDF coordinate frame/order mismatch')
            origin, origin_method, origin_sources = origin_for(folder, arrays)
            target = arrays['initial_coords_A'] + origin
            fixed = task['fixed_sampler_indices']
            require(np.max(np.abs(target[fixed] - ref.GetConformer().GetPositions()[task['fixed_reference_indices']])) < .001,
                    'Prior/reference fixed coordinate discrepancy')
            protein_hashes = {h for p, h in info['input_sha256'].items() if Path(p).suffix.lower() == '.pdb'}
            candidates = [root / 'inputs/protein.pdb', group / 'inputs/protein.pdb']
            candidates += [Path(p) for p in info['input_sha256'] if Path(p).suffix.lower() == '.pdb']
            if protein_override:
                candidates.insert(0, Path(protein_override))
            proteins = [p for p in candidates if p.is_file() and sha(p) in protein_hashes]
            require(proteins, 'Original full protein PDB missing. Supply the same file with PROTEIN_PDB; do not substitute a pocket.')
            protein = proteins[0]
            trajectory = {}
            for stream, metadata in report['streams'].items():
                trajectory[stream + '_coords'] = arrays[stream + '_coords'].astype(float) + origin.astype(float)
                trajectory[stream + '_atom_classes'] = arrays[stream + '_atom_classes']
                trajectory[stream + '_valid'] = arrays[stream + '_valid']
            require(np.array_equal(arrays['sampler_state_atom_classes'][-1][valid], classes[valid]), 'Terminal atom identities differ')
            files = [quality_path, sdf, reference_path, task_path, array_path, report_path, protein, *origin_sources]
            if study_path.is_file():
                files.append(study_path)
            items.append(dict(id=rid, input_id=f"{info['case']}_lig{info.get('ligand_idx', 'unknown')}",
                seed=info['seed'], mol=mol, ref=ref, protein=protein, raw_sdf=sdf.read_bytes(),
                xyz=xyz, target=target, slots=slots, task=task, expected=rows[rid],
                raw_fixed_evidence=info.get('raw_graph'),
                trajectory=trajectory, trajectory_meta=dict(streams=report['streams'], atom_vocabulary=vocab),
                imported_observations=dict(connection_report={k: report.get(k) for k in
                    ('diagnosis', 'decoded_checks', 'timing_note', 'scope')}, study=study.get(rid),
                    status='preserved_prior_observations_not_new_neural_analysis'),
                provenance=dict(adapter='flowr_records', source_root=str(root), origin_method=origin_method,
                    source_hashes={str(p): sha(p) for p in files}, source_versions=quality['versions'],
                    model_identity=info.get('recorded_source_sha256'), checkpoint_sha256=info.get('checkpoint_sha256'),
                    fixed_definition='complement_of_reference_Murcko_scaffold')))
    require(items and len({r['id'] for r in items}) == len(items), 'No samples or duplicate sample IDs')
    return items

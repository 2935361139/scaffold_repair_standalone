"""Model-neutral SDF/manifest input with explicit, zero-based atom mappings."""
from pathlib import Path
import numpy as np
from rdkit import Chem
from ..io import read_json, require, resolve, sha, safe_id
from ..evaluation import task_from_mapping


def load_one(spec, base):
    sdf = resolve(spec['sdf'], base)
    ref_path = resolve(spec['reference_sdf'], base)
    protein = resolve(spec['protein_pdb'], base)
    for path in (sdf, ref_path, protein):
        require(path.is_file(), f'Input file missing: {path}')
    supplier = Chem.SDMolSupplier(str(sdf), sanitize=False, removeHs=False)
    index = spec.get('sdf_index', 0)
    require(type(index) is int and 0 <= index < len(supplier), 'SDF record index out of range')
    mol = supplier[index]
    refs = Chem.SDMolSupplier(str(ref_path), sanitize=False, removeHs=False)
    reference_index = spec.get('reference_index', 0)
    require(type(reference_index) is int and 0 <= reference_index < len(refs), 'Reference record index out of range')
    ref = refs[reference_index]
    require(ref is not None, 'Reference molecule unreadable')
    Chem.SanitizeMol(ref)
    require(ref.GetNumConformers() == 1 and np.isfinite(ref.GetConformer().GetPositions()).all(), 'Reference needs finite coordinates')
    require(all(a.GetAtomicNum() > 1 for a in ref.GetAtoms()), 'Use heavy-atom reference SDF; do not silently renumber atom mappings')
    fixed = spec['fixed_generated_indices']
    ref_fixed = spec['fixed_reference_indices']
    task = task_from_mapping(ref, fixed, ref_fixed)
    n = mol.GetNumAtoms() if mol is not None else 0
    if mol is not None:
        require(mol.GetNumConformers() == 1, 'Generated molecule needs one conformer')
        require(all(a.GetAtomicNum() > 1 for a in mol.GetAtoms()), 'Use heavy-atom generated SDF and its matching atom indices')
        require(all(i < n for i in fixed), 'Fixed generated atom index out of range')
        xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=np.float64)
    else:
        xyz = np.zeros((0, 3), dtype=np.float64)
    sources = [sdf, ref_path, protein]
    if spec.get('coordinates_npy'):
        path = resolve(spec['coordinates_npy'], base)
        accurate = np.load(path, allow_pickle=False)
        require(accurate.shape == (n, 3) and accurate.dtype in (np.dtype('float32'), np.dtype('float64')),
                'High precision coordinate array must be float32/float64, shape (atoms, 3)')
        require(np.isfinite(accurate).all() and (not n or np.max(np.abs(accurate - xyz)) <= .002),
                'High precision coordinates differ from SDF atom order/frame')
        xyz = accurate
        sources.append(path)
    require(np.isfinite(xyz).all(), 'Nonfinite generated coordinates')
    target = xyz.copy()
    if mol is not None:
        target[fixed] = ref.GetConformer().GetPositions()[ref_fixed]
    if spec.get('target_coordinates_npy'):
        path = resolve(spec['target_coordinates_npy'], base)
        target = np.load(path, allow_pickle=False)
        require(target.shape == xyz.shape and target.dtype == xyz.dtype and np.isfinite(target).all(),
                'Target coordinate shape/dtype mismatch')
        require(not n or np.max(np.abs(target[fixed] - ref.GetConformer().GetPositions()[ref_fixed])) < .001,
                'Target fixed coordinates differ from reference')
        sources.append(path)
    trajectory, metadata = {}, {}
    require(bool(spec.get('trajectory_npz')) == bool(spec.get('trajectory_metadata')), 'Trajectory needs both NPZ and metadata JSON')
    if spec.get('trajectory_npz'):
        path, meta_path = resolve(spec['trajectory_npz'], base), resolve(spec['trajectory_metadata'], base)
        with np.load(path, allow_pickle=False) as arrays:
            trajectory = {k: arrays[k] for k in arrays.files}
        metadata = read_json(meta_path)
        sources += [path, meta_path]
        require(metadata.get('coordinate_frame') == 'protein_world_angstrom', 'Trajectory must already be in protein coordinates, in angstroms')
        require(mol is not None, 'Cannot map a trajectory to an unreadable endpoint')
        for stream, stream_meta in metadata['streams'].items():
            coords = trajectory[stream + '_coords']
            require(coords.ndim == 3 and coords.shape[1:] == (n, 3) and np.isfinite(coords).all(), 'Trajectory atom count/frame error')
            require(len(stream_meta['indices']) == len(coords), 'Trajectory metadata length mismatch')
        terminal = metadata['terminal_stream']
        require(np.allclose(trajectory[terminal + '_coords'][-1], xyz, atol=2e-5, rtol=0), 'Trajectory endpoint does not match molecule')
        require(all(metadata['atom_vocabulary'][int(c)] == a.GetSymbol() for c, a in
                    zip(trajectory[terminal + '_atom_classes'][-1], mol.GetAtoms())), 'Trajectory terminal element order mismatch')
        require(bool(np.all(trajectory[terminal + '_valid'][-1])), 'Trajectory endpoint masks retained atoms')
    return dict(id=safe_id(spec['id']), input_id=str(spec['input_id']), seed=spec.get('seed'),
        mol=mol, ref=ref, protein=protein, raw_sdf=supplier.GetItemText(index).encode('utf-8'),
        xyz=xyz, target=target, slots=list(range(n)), task=task, expected=None, raw_fixed_evidence=None,
        trajectory=trajectory, trajectory_meta=metadata, imported_observations=None,
        provenance=dict(adapter='sdf', model=str(spec.get('model', 'external')), sdf_index=index,
            reference_index=reference_index, source_hashes={str(p): sha(p) for p in sources},
            fixed_definition='explicit_user_mapping',
            precision='original_array' if spec.get('coordinates_npy') else 'SDF_coordinates_only',
            warning='SDF rounding may change optimization near thresholds; exact prior equivalence requires original coordinate arrays.'))


def load_manifest(path):
    path = Path(path).resolve()
    data = read_json(path)
    require(data.get('schema') == 'scaffold_repair_manifest_v1', 'Unknown manifest schema')
    require(isinstance(data.get('samples'), list) and data['samples'], 'Manifest contains no samples')
    items = [load_one(spec, path.parent) for spec in data['samples']]
    for item in items:
        item['provenance']['source_hashes'][str(path)] = sha(path)
    require(len({r['id'] for r in items}) == len(items), 'Duplicate sample IDs')
    return items


def load_sdf(sdf_path, protein_path, reference_path, mapping_path, record_indices=None):
    path = Path(mapping_path).resolve()
    data = read_json(path)
    require(data.get('schema') == 'scaffold_repair_mapping_v1', 'Unknown atom mapping schema')
    supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    indices = list(range(len(supplier))) if record_indices is None else record_indices
    require(indices and len(set(indices)) == len(indices), 'Empty or duplicated SDF record selection')
    specs = data['samples']
    require(len({s['sdf_index'] for s in specs}) == len(specs), 'Duplicate SDF indices in mapping file')
    by_index = {s['sdf_index']: s for s in specs}
    require(set(indices) <= set(by_index), 'Every selected SDF record needs its own explicit atom mapping')
    items = []
    for index in indices:
        spec = dict(by_index[index], sdf=str(sdf_path), protein_pdb=str(protein_path), reference_sdf=str(reference_path))
        item = load_one(spec, path.parent)
        item['provenance']['source_hashes'][str(path)] = sha(path)
        items.append(item)
    require(len({r['id'] for r in items}) == len(items), 'Duplicate sample IDs')
    return items

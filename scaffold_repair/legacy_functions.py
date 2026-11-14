"""Unchanged scientific functions extracted from verified previous sources."""
import numpy as np
from rdkit import Chem
from .geometry_intervention import require

def coordinate_arm(coords, target, fixed, kind, seed):
    import numpy as np
    x=np.asarray(coords); goal=np.asarray(target)
    require(x.shape==goal.shape and x.ndim==2 and x.shape[1]==3, 'Invalid coordinate shape')
    require(x.dtype.kind=='f' and np.isfinite(x).all() and np.isfinite(goal).all(), 'Invalid coordinates')
    require(len(fixed)>0 and len(set(fixed))==len(fixed) and all(0<=i<len(x) for i in fixed), 'Invalid fixed selection')
    require(kind in ('restore_fixed','random_direction'), 'Unknown coordinate arm')
    result=x.copy(); delta=goal[fixed].astype(float)-x[fixed].astype(float)
    if kind=='restore_fixed':
        result[fixed]=goal[fixed]
    else:
        rng=np.random.default_rng(seed); change=np.zeros_like(delta)
        for k,d in enumerate(delta):
            length=float(np.linalg.norm(d))
            if length==0:continue
            direction=d/length
            for _ in range(16):
                noise=rng.normal(size=3); noise-=np.dot(noise,direction)*direction
                norm=float(np.linalg.norm(noise))
                if norm>1e-12:break
            else:raise ValueError('Degenerate random direction')
            change[k]=noise*(length/norm)
        result[fixed]=(x[fixed].astype(float)+change).astype(x.dtype)
    actual=result[fixed].astype(float)-x[fixed].astype(float)
    norms=np.linalg.norm(delta,axis=1);actual_norms=np.linalg.norm(actual,axis=1)
    tolerance=16*np.finfo(x.dtype).eps*(np.linalg.norm(x[fixed].astype(float),axis=1)+norms+1.)
    require(np.all(np.abs(norms-actual_norms)<=tolerance), 'Control magnitude mismatch')
    other=np.ones(len(x),dtype=bool);other[fixed]=False
    require(np.array_equal(result[other],x[other]), 'Unselected coordinates changed')
    require(np.isfinite(result).all(), 'Nonfinite patched coordinates')
    informative=norms>tolerance
    cosines=np.divide(np.sum(actual*delta,axis=1),actual_norms*norms,
                      out=np.zeros_like(norms),where=(actual_norms*norms)>0)
    if kind=='random_direction':
        require(np.all(np.abs(cosines[informative])<=tolerance[informative]/actual_norms[informative]),
                'Control direction not orthogonal within numerical precision')
    return result,dict(kind=kind,seed=seed if kind=='random_direction' else None,
        cosines_to_target=cosines.tolist(),
        selected_slots=list(fixed),target_norms_A=norms.tolist(),actual_norms_A=actual_norms.tolist(),
        norm_rounding_tolerance_A=tolerance.tolist(),informative=(norms>tolerance).tolist(),
        only_selected_coordinates_changed=True)

def fixed_diagnostics(mol, ref, task, slots, np):
    fixed = list(task["fixed_sampler_indices"])
    ref_map = dict(zip(fixed, task["fixed_reference_indices"]))
    gen_map = {slot: index for index, slot in enumerate(slots)}
    if any(i not in gen_map for i in fixed):
        raise ValueError("A recorded fixed atom is absent; requires a separate mapping audit.")
    expected = {i: [] for i in fixed}
    for entry in task["reference_attachment_sites"]:
        expected[entry["fixed_sampler_index"]].append(entry["reference_bond_type"])
    actual = {i: [] for i in fixed}
    neighbors = {i: [] for i in fixed}
    for bond in mol.GetBonds():
        a, b = slots[bond.GetBeginAtomIdx()], slots[bond.GetEndAtomIdx()]
        for i, j in [(a, b), (b, a)]:
            if i in actual and j not in actual:
                actual[i].append(str(bond.GetBondType()))
                neighbors[i].append(j)
    violations = [i for i in fixed if sorted(actual[i]) != sorted(expected[i])]
    elements = [i for i in fixed if mol.GetAtomWithIdx(gen_map[i]).GetAtomicNum()
                != ref.GetAtomWithIdx(ref_map[i]).GetAtomicNum()]
    charges = [i for i in fixed if mol.GetAtomWithIdx(gen_map[i]).GetFormalCharge()
               != ref.GetAtomWithIdx(ref_map[i]).GetFormalCharge()]
    def bond_type(m, i, j):
        bond = m.GetBondBetweenAtoms(i, j)
        return str(bond.GetBondType()) if bond else None
    internal = [[i, j] for k, i in enumerate(fixed) for j in fixed[:k]
                if bond_type(mol, gen_map[i], gen_map[j]) != bond_type(ref, ref_map[i], ref_map[j])]
    delta = np.asarray([mol.GetConformer().GetPositions()[gen_map[i]]
                        - ref.GetConformer().GetPositions()[ref_map[i]] for i in fixed])
    distances = np.linalg.norm(delta, axis=1)
    return {
        "fixed_element_changed_sites": elements,
        "fixed_charge_changed_sites": charges,
        "fixed_internal_bond_changed_pairs": internal,
        "fixed_rms_displacement_A": float(np.sqrt(np.mean(distances**2))),
        "fixed_max_displacement_A": float(distances.max()),
        "strict_boundary_reference_match": not violations,
        "strict_boundary_differing_sites": violations,
        "attachment_sites_without_generated_neighbor": [i for i in fixed if expected[i] and not actual[i]],
        "nonattachment_sites_with_generated_neighbors": [i for i in fixed if not expected[i] and actual[i]],
        "boundary_details": [{"fixed_site": i, "reference_bond_types": sorted(expected[i]),
                              "generated_bond_types": sorted(actual[i]), "generated_neighbors": neighbors[i]}
                             for i in fixed],
    }

def pb_metrics(mol, protein, PoseBusters, Chem):
    # Matches generate_from_pdb's unoptimized heavy-atom model path: sanitize,
    # add H coordinates, then evaluate_pb_validity's PoseBusters(config="dock").
    import numpy as np
    with_hs = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    buster = PoseBusters(config="dock")
    table = buster.bust([with_hs], None, str(protein))
    if len(table) != 1 or not len(table.columns):
        raise RuntimeError("PoseBusters returned no unique, nonempty result row.")
    row = table.iloc[0]
    if row.isna().any():
        raise RuntimeError("PoseBusters returned missing check results; not a molecular failure label.")
    if any(not isinstance(v, (bool, np.bool_)) for v in row.tolist()):
        raise RuntimeError("Unexpected nonboolean PoseBusters result columns.")
    checks = {str(k): bool(v) for k, v in row.items()}
    return {"posebusters_valid": bool(table.all(axis=1).iloc[0]),
            "posebusters_checks": checks,
            "posebusters_failed_checks": [k for k, v in checks.items() if not v]}

def _sanitize_or_fix_aromaticity(mol: Chem.Mol) -> bool:
    """Try to sanitize *mol* in-place.

    When molecules are reconstructed via ``to_rdkit(sanitise=False)``, bonds
    may carry the ``AROMATIC`` type while their atoms lack the corresponding
    aromatic flag.  ``SanitizeMol`` then fails during kekulization because it
    cannot resolve the aromatic system.

    If the first sanitization attempt fails, this helper propagates aromaticity
    flags from bonds to atoms and retries.  Returns ``True`` on success.
    """
    try:
        Chem.SanitizeMol(mol)
        return True
    except Exception:
        pass

    # Fix aromatic inconsistency: bonds may have AROMATIC type but atoms lack
    # the aromatic flag.  Propagate the flag so kekulization can succeed.
    try:
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.GetBeginAtom().SetIsAromatic(True)
                bond.GetEndAtom().SetIsAromatic(True)
                bond.SetIsAromatic(True)
        Chem.SanitizeMol(mol)
        return True
    except Exception:
        return False

def suffix(points):
    require(bool(points), 'Empty observation series')
    first = next((p for p in points if p['value'] is True), None)
    if points[-1]['value'] is not True:
        return dict(status='terminal_signature_absent_or_unknown', first_observed=first,
                    first_persistent_observation=None, causal_onset=None)
    i = len(points)-1
    while i and points[i-1]['value'] is True:
        i -= 1
    return dict(status='observed', first_observed=first,
                first_persistent_observation=points[i], last_observation=points[-1],
                preceding_observation=points[i-1] if i else None,
                boundary_is_confirmed_change=i > 0 and points[i-1]['value'] is False,
                causal_onset=None)

def series(xyz, classes, valid, vocab, slot, element, environment_xyz,
           radius_sum_scaled, cutoff, metadata):
    n = len(xyz)
    require(metadata['indices'] == list(range(n)), 'Nonconsecutive temporal indices')
    require(len(metadata['times']) == len(metadata['stages']) == n, 'Incomplete time metadata')
    require(classes.shape == valid.shape == xyz.shape[:2] and xyz.shape[-1] == 3,
            'Trajectory array shape mismatch')
    require(0 <= slot < xyz.shape[1], 'Unmapped clash slot')
    require(np.isfinite(xyz).all() and radius_sum_scaled > 0, 'Invalid coordinates/radii')
    points = []
    for k in range(n):
        cls = int(classes[k, slot])
        require(0 <= cls < len(vocab), 'Invalid atom vocabulary index')
        same = bool(valid[k, slot]) and vocab[cls] == element
        distance = float(np.linalg.norm(xyz[k, slot] - environment_xyz))
        # Numerical ambiguity is unknown, never a new molecular failure threshold.
        margin = distance - radius_sum_scaled * cutoff
        value = bool(margin < 0) if same and abs(margin) > 1e-5 else None
        points.append(dict(index=k, native_time=metadata['times'][k], stage=metadata['stages'][k],
                           value=value, distance_A=distance,
                           relative_distance=distance/radius_sum_scaled,
                           endpoint_element_present=same))
    return points


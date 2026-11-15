"""Prepare a receptor copy, removing only an explicitly identified reference ligand.

No alignment, arbitrary clash deletion, water removal, or endpoint-dependent selection.
"""
from collections import defaultdict
from pathlib import Path
import numpy as np
from rdkit import Chem
from .io import require, sha, write_json

DEFAULTS = dict(auto_exclude_reference=True, exclude_residues=[], match_tolerance_A=0.01)


def residue_key(info):
    return (info.GetChainId().strip(), info.GetResidueNumber(), info.GetInsertionCode().strip(),
            info.GetResidueName().strip())


def prepare(source, reference, destination, options=None):
    options = {**DEFAULTS, **(options or {})}
    require(type(options['auto_exclude_reference']) is bool, 'auto_exclude_reference must be boolean')
    tolerance = options['match_tolerance_A']
    require(isinstance(tolerance, (float, int)) and 0 < tolerance <= .1, 'Invalid reference matching tolerance')
    source, destination = Path(source), Path(destination)
    text = source.read_text()
    require(sum(line.startswith('MODEL ') for line in text.splitlines()) <= 1,
            'Multi-model receptor PDB is ambiguous; supply one model')
    mol = Chem.MolFromPDBBlock(text, sanitize=False, removeHs=False, proximityBonding=False)
    require(mol is not None and mol.GetNumConformers() == 1, 'Unreadable receptor')
    groups = defaultdict(list)
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        require(info is not None, 'Receptor needs PDB residue identifiers')
        groups[residue_key(info)].append(atom.GetIdx())
    explicit = set()
    for spec in options['exclude_residues']:
        require(isinstance(spec, dict) and {'chain', 'number', 'name'} <= spec.keys(),
                'Each excluded residue needs chain, number, name, optional insertion')
        key = (str(spec['chain']), int(spec['number']), str(spec.get('insertion', '')), str(spec['name']))
        require(key in groups, f'Excluded residue not found: {key}')
        require(all(mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetIsHeteroAtom() for i in groups[key]),
                'Cannot exclude protein ATOM residues as a repair')
        explicit.add(key)
    matches = []
    target = reference.GetConformer().GetPositions()
    elems = [a.GetAtomicNum() for a in reference.GetAtoms()]
    require(all(e > 1 for e in elems), 'Reference receptor matching requires heavy atoms')
    for key, ids in groups.items():
        if key[3] in {'HOH', 'WAT', 'DOD'} or not all(mol.GetAtomWithIdx(i).GetPDBResidueInfo().GetIsHeteroAtom() for i in ids):
            continue
        heavy = [i for i in ids if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
        if len(heavy) != len(elems):
            continue
        xyz = mol.GetConformer().GetPositions()[heavy]
        d = np.linalg.norm(target[:, None, :] - xyz[None, :, :], axis=-1)
        valid = (d <= tolerance) & (np.array(elems)[:, None] == np.array([mol.GetAtomWithIdx(i).GetAtomicNum() for i in heavy])[None, :])
        if np.all(valid.sum(axis=0) == 1) and np.all(valid.sum(axis=1) == 1):
            matches.append(key)
    automatic = set()
    if options['auto_exclude_reference']:
        remaining = [k for k in matches if k not in explicit]
        require(len(remaining) <= 1, 'Multiple reference ligand residues match; explicitly select exclusions')
        automatic.update(remaining)
    removed = automatic | explicit
    # Preserve all retained PDB coordinate records verbatim; remove dangling CONECT references.
    kept, serials, removed_lines = [], set(), 0
    for line in text.splitlines(keepends=True):
        if line.startswith(('ATOM  ', 'HETATM')):
            key = (line[21:22].strip(), int(line[22:26]), line[26:27].strip(), line[17:20].strip())
            if key in removed:
                removed_lines += 1
                continue
            serials.add(int(line[6:11]))
        if not line.startswith('CONECT'):
            kept.append(line)
    # Protein loading uses proximityBonding=False. Retain connectivity between retained atoms.
    connects = []
    for line in text.splitlines():
        if line.startswith('CONECT'):
            nums = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            if nums and nums[0] in serials:
                neighbors = [n for n in nums[1:] if n in serials]
                if neighbors:
                    connects.append('CONECT' + ''.join(f'{n:5d}' for n in [nums[0], *neighbors]) + '\n')
    # Place connectivity before END so loaders do not ignore it.
    kept = [line for line in kept if line.strip() != 'END']
    if not removed:
        destination.write_bytes(source.read_bytes())
    else:
        destination.write_text(''.join(kept + connects) + 'END\n')
    result = Chem.MolFromPDBFile(str(destination), sanitize=False, removeHs=False, proximityBonding=False)
    require(result is not None and result.GetNumAtoms() > 0, 'Empty prepared receptor')
    require(any(not a.GetPDBResidueInfo().GetIsHeteroAtom() for a in result.GetAtoms()), 'Receptor contains no protein ATOM records')
    audit = dict(source_sha256=sha(source), prepared_sha256=sha(destination), settings=options,
                 matched_reference_residues=[list(k) for k in matches],
                 automatic_exclusions=[list(k) for k in sorted(automatic)],
                 explicit_exclusions=[list(k) for k in sorted(explicit)], removed_coordinate_records=removed_lines,
                 unchanged_retained_coordinates=True, is_molecular_repair=False)
    write_json(destination.with_name('receptor_preparation.json'), audit)
    return audit

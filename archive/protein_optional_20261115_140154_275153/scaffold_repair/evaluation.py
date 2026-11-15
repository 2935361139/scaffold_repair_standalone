"""Original endpoint conjunction; explicit fixed-atom mapping for any model."""
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from posebusters import PoseBusters
from .io import require
from .import_state import restore as restore_import_state, encode as encode_import_state, PROPERTY as IMPORT_STATE_PROPERTY
from .legacy_functions import fixed_diagnostics, pb_metrics, _sanitize_or_fix_aromaticity

ENDPOINT_KEYS = ('rdkit_valid', 'fc_valid', 'root_fixed_substructure_match',
                 'fixed_chemistry_preserved', 'posebusters_valid')
COMPARISON_KEYS = (*ENDPOINT_KEYS, 'posebusters_checks', 'murcko_scaffold_changed',
                   'molecule_smiles', 'scaffold_smiles')


def read_mol(path):
    if not path.stat().st_size:
        return None
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    require(len(supplier) == 1, f'Expected one SDF record: {path}')
    return restore_import_state(supplier[0])


def restore_flowr_decoder_representation(mol):
    """The original tensor decoder did not assign chirality from 3D coordinates.

    RDKit's SDF reader does so even with sanitize=False. Undo that reader-only
    assignment for FLOWR records; external SDF stereochemistry is left intact.
    """
    if mol is not None:
        Chem.RemoveStereochemistry(mol)
        if mol.HasProp('recorded_sample_index'):
            mol.ClearProp('recorded_sample_index')
    return mol


def with_coords(mol, coords):
    require(mol is not None and mol.GetNumConformers() == 1, 'Expected a molecule with one conformer')
    require(coords.shape == (mol.GetNumAtoms(), 3) and np.isfinite(coords).all(), 'Invalid coordinates')
    result = Chem.Mol(mol)
    for i, xyz in enumerate(coords):
        result.GetConformer().SetAtomPosition(i, tuple(float(v) for v in xyz))
    return result


def save_mol(path, mol):
    with Chem.SDWriter(str(path)) as writer:
        if mol is not None:
            if mol.HasProp(IMPORT_STATE_PROPERTY):
                writer.SetKekulize(False)
                writer.SetForceV3000(True)
                encode_import_state(mol)
            writer.write(mol)


def expected_fixed(ref):
    molecule = Chem.Mol(ref)
    require(_sanitize_or_fix_aromaticity(molecule), 'Reference molecule cannot be sanitized')
    for atom in molecule.GetAtoms():
        atom.SetIntProp('org_idx', atom.GetIdx())
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    removed = {a.GetIntProp('org_idx') for a in scaffold.GetAtoms()}
    return [i for i in range(molecule.GetNumAtoms()) if i not in removed]


def fixed_substructure_match(mol, ref, indices):
    reference = Chem.Mol(ref)
    require(_sanitize_or_fix_aromaticity(reference), 'Invalid reference chemistry')
    query = Chem.RWMol(reference)
    for i in reversed(range(reference.GetNumAtoms())):
        if i not in indices:
            query.RemoveAtom(i)
    return bool(mol.HasSubstructMatch(query.GetMol()))


def task_from_mapping(ref, fixed, ref_fixed):
    require(len(fixed) == len(ref_fixed) > 0, 'Fixed mapping must contain at least one atom pair')
    require(len(set(fixed)) == len(fixed) and len(set(ref_fixed)) == len(ref_fixed), 'Duplicate fixed atom mapping')
    require(all(type(i) is int and i >= 0 for i in fixed), 'Invalid generated atom index')
    require(all(type(i) is int and 0 <= i < ref.GetNumAtoms() for i in ref_fixed), 'Invalid reference atom index')
    mapping = dict(zip(ref_fixed, fixed))
    attachments = []
    for bond in ref.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for i, j in ((a, b), (b, a)):
            if i in mapping and j not in mapping:
                attachments.append(dict(fixed_reference_index=i, fixed_sampler_index=mapping[i],
                                        removed_scaffold_reference_index=j, reference_bond_type=str(bond.GetBondType())))
    return dict(fixed_sampler_indices=fixed, fixed_reference_indices=ref_fixed,
                reference_attachment_sites=attachments)


def evaluate(mol, coords, ref, task, slots, protein, raw_fixed_evidence=None):
    row = dict.fromkeys(ENDPOINT_KEYS)
    row.update(rdkit_valid=False, fc_valid=False, posebusters_checks=None,
               posebusters_failed_checks=[], murcko_scaffold_changed=None,
               molecule_smiles=None, scaffold_smiles=None)
    row['molecule_readable'] = mol is not None
    if raw_fixed_evidence is not None:
        row['fixed_chemistry_preserved'] = raw_fixed_evidence.get('fixed_chemistry_preserved')
        row['raw_fixed_evidence_source'] = 'previous_checksum_verified_raw_graph_audit'
    if mol is None:
        row['geometry_status'] = 'not_assessed_unreadable_generated_molecule'
        return row, None
    result = with_coords(mol, coords)
    row['fixed_missing_sites'] = sorted(set(task['fixed_sampler_indices']) - set(slots))
    if not row['fixed_missing_sites']:
        # Fixed identity / bond discrepancies can be reported even before the
        # complete graph passes RDKit sanitization.
        row.update(fixed_diagnostics(result, ref, task, slots, np))
    try:
        Chem.SanitizeMol(result)
    except Exception:
        row['geometry_status'] = 'not_assessed_invalid_chemistry'
        return row, result
    row['rdkit_valid'] = True
    row['fc_valid'] = len(Chem.GetMolFrags(result)) == 1
    fixed = task['fixed_sampler_indices']
    missing = sorted(set(fixed) - set(slots))
    row['fixed_missing_sites'] = missing
    if missing:
        row['fixed_chemistry_preserved'] = False
    else:
        row.update(fixed_diagnostics(result, ref, task, slots, np))
        row['fixed_chemistry_preserved'] = not any(row[key] for key in
            ('fixed_element_changed_sites', 'fixed_charge_changed_sites', 'fixed_internal_bond_changed_pairs'))
    row.pop('raw_fixed_evidence_source', None)
    row['root_fixed_substructure_match'] = fixed_substructure_match(result, ref, task['fixed_reference_indices'])
    row['molecule_smiles'] = Chem.MolToSmiles(result, isomericSmiles=False)
    row['scaffold_smiles'] = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(result), isomericSmiles=False)
    row['murcko_scaffold_changed'] = row['scaffold_smiles'] != Chem.MolToSmiles(
        MurckoScaffold.GetScaffoldForMol(ref), isomericSmiles=False)
    row.update(pb_metrics(result, protein, PoseBusters, Chem))
    row['geometry_status'] = 'evaluated'
    return row, result


def compare_endpoint(actual, expected):
    for key in COMPARISON_KEYS:
        require(actual.get(key) == expected.get(key),
                f'Baseline mismatch in {key}: new={actual.get(key)!r}, previous={expected.get(key)!r}')
    for key in ('fixed_rms_displacement_A', 'fixed_max_displacement_A'):
        if key in expected:
            require(key in actual and abs(actual[key] - expected[key]) <= 1e-5,
                    f'Baseline displacement mismatch in {key}')

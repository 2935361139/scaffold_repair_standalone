"""Protein related functions."""
from __future__ import annotations

import logging
from typing import Iterable

from rdkit.Chem.rdchem import Atom, Mol

logger = logging.getLogger(__name__)
_inorganic_cofactor_elements = {
    "Li",
    "Be",
    "Na",
    "Mg",
    "Cl",
    "K",
    "Ca",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Br",
    "Rb",
    "Mo",
    "Cd",
}
_inorganic_cofactor_ccd_codes = {
    "FES",
    "MOS",
    "PO3",
    "PO4",
    "PPK",
    "SO3",
    "SO4",
    "VO4",
}
_water_ccd_codes = {
    "HOH",
}


def get_atom_type_mask(mol: Mol, ignore_types: Iterable[str]) -> list[bool]:
    """Get mask for atoms to keep."""
    ignore_types = set(ignore_types)
    unsupported = ignore_types - {"hydrogens", "protein", "organic_cofactors", "inorganic_cofactors", "waters"}
    if unsupported:
        raise ValueError(f"Ignore types {unsupported} not supported.")

    if not mol.GetNumAtoms():
        return []
    if mol.GetAtomWithIdx(0).GetPDBResidueInfo() is None:
        logger.warning("No PDB information found. Assuming organic molecule.")
    sizes = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and atom.GetAtomicNum() > 1:
            key = (info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode(), info.GetResidueName())
            sizes[key] = sizes.get(key, 0) + 1
    result = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            category = 'organic_cofactors'
        elif not info.GetIsHeteroAtom():
            category = 'protein'
        elif info.GetResidueName().strip() in _water_ccd_codes | {'WAT', 'DOD'}:
            category = 'waters'
        else:
            key = (info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode(), info.GetResidueName())
            inorganic = info.GetResidueName().strip() in _inorganic_cofactor_ccd_codes or (
                sizes.get(key, 0) == 1 and atom.GetSymbol() in _inorganic_cofactor_elements)
            category = 'inorganic_cofactors' if inorganic else 'organic_cofactors'
        result.append(not (('hydrogens' in ignore_types and atom.GetAtomicNum() == 1) or category in ignore_types))
    return result


def _keep_atom(  # noqa: PLR0913, PLR0911
    atom: Atom, ignore_h: bool, ignore_protein: bool, ignore_org_cof: bool, ignore_inorg_cof: bool, ignore_water: bool
) -> bool:
    """Whether to keep atom for given ignore flags."""
    symbol = atom.GetSymbol()
    if ignore_h and symbol == "H":
        return False

    if ignore_inorg_cof and symbol in _inorganic_cofactor_elements:
        return False

    # if loaded from PDB file, we can use the residue names and the hetero flag
    info = atom.GetPDBResidueInfo()
    if info is None:
        if ignore_org_cof:
            return False
        return True

    is_hetero = info.GetIsHeteroAtom()
    if ignore_protein and not is_hetero:
        return False

    residue_name = info.GetResidueName()
    if ignore_water and residue_name in _water_ccd_codes:
        return False

    if ignore_inorg_cof and residue_name in _inorganic_cofactor_ccd_codes:
        return False

    return True

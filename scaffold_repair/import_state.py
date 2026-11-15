"""Preserve structural atom state across SDF serialization; never accept labels."""
import json
from rdkit import Chem

PROPERTY = 'V2_IMPORT_ATOM_STATE_V1'


def bond_signature(mol):
    # SDF writers may reverse the endpoints of an undirected covalent bond
    # when encoding stereochemistry. The molecular graph is unchanged.
    return sorted([[min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                    max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), str(b.GetBondType())]
                   for b in mol.GetBonds()])


def encode(mol):
    # Explicit H count matters for invalid valence and aromatic [nH].
    data = {'atoms': [[a.GetAtomicNum(), a.GetFormalCharge(), a.GetNumExplicitHs(),
                       a.GetNoImplicit(), a.GetIsAromatic()] for a in mol.GetAtoms()],
            'bonds': bond_signature(mol)}
    mol.SetProp(PROPERTY, json.dumps(data, separators=(',', ':')))
    return mol


def restore(mol):
    if mol is None or not mol.HasProp(PROPERTY):
        return mol
    data = json.loads(mol.GetProp(PROPERTY))
    expected_bonds = sorted([[min(i,j), max(i,j), kind] for i,j,kind in data['bonds']])
    if len(data['atoms']) != mol.GetNumAtoms() or expected_bonds != bond_signature(mol):
        raise ValueError('Imported atom state does not match SDF graph')
    for atom, row in zip(mol.GetAtoms(), data['atoms']):
        z, charge, nh, no_implicit, aromatic = row
        if z != atom.GetAtomicNum() or charge != atom.GetFormalCharge() or type(nh)is not int or nh<0 or nh>64 or type(no_implicit)is not bool or type(aromatic)is not bool:
            raise ValueError('Invalid structural atom state')
        atom.SetNumExplicitHs(nh)
        atom.SetNoImplicit(no_implicit)
        atom.SetIsAromatic(aromatic)
    mol.UpdatePropertyCache(strict=False)
    return mol


def block(mol):
    mol = encode(Chem.Mol(mol))
    # A real SDF writer retains properties; no sanitize/kekulize step.
    import io
    stream=io.StringIO()
    writer=Chem.SDWriter(stream)
    writer.SetKekulize(False)
    writer.SetForceV3000(True)
    writer.write(mol)
    writer.close()
    return stream.getvalue()

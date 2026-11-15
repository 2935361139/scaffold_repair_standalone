"""Local hypotheses for RDKit-invalid graphs, without atom/charge replacement."""
from rdkit import Chem, rdBase


def problems(mol):
    if mol is None:
        return [dict(type='UnreadableMolecule', atoms=[], message='No molecular graph was loaded')]
    try:
        with rdBase.BlockLogs():
            found = Chem.DetectChemistryProblems(Chem.Mol(mol))
        return [dict(type=p.GetType(), atoms=([p.GetAtomIdx()] if hasattr(p,'GetAtomIdx') else
            list(p.GetAtomIndices()) if hasattr(p,'GetAtomIndices') else []), message=p.Message()) for p in found]
    except Exception as exc:
        return [dict(type='ChemistryInspectionError', atoms=[], message=str(exc))]


def repair_invalid(mol, ref, mapping, budget, max_edits=2):
    """Reduce excessive bond order / remove inconsistent edges around flagged atoms.

    Fixed-induced bonds and all formal charges/elements are immutable here.
    This does not guess an intended generated-core bond order when it is valid.
    """
    from .topology import graph_key, fixed_matches, LIMITS
    if not problems(mol):
        return []
    pending, seen, complete = [Chem.Mol(mol)], {graph_key(mol)}, []
    for _ in range(max_edits):
        next_states = []
        for base in pending:
            bad = {i for p in problems(base) for i in p['atoms']}
            for bond in sorted(base.GetBonds(), key=lambda b:tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx())))):
                a,b=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
                if not ({a,b} & bad) or (a in mapping and b in mapping):
                    continue
                orders = {Chem.BondType.TRIPLE:[Chem.BondType.DOUBLE,Chem.BondType.SINGLE],
                          Chem.BondType.DOUBLE:[Chem.BondType.SINGLE],
                          Chem.BondType.SINGLE:[None]}.get(bond.GetBondType(),[])
                for order in orders:
                    if budget.attempts >= LIMITS['max_sanitizations']:
                        return complete
                    changed=Chem.RWMol(base)
                    if order is None: changed.RemoveBond(a,b)
                    else: changed.GetBondBetweenAtoms(a,b).SetBondType(order)
                    key=graph_key(changed)
                    if key in seen: continue
                    seen.add(key)
                    valid=budget.sanitize(changed)
                    if valid is not None and fixed_matches(valid,ref,mapping):
                        complete.append(valid)
                        if len(complete)>=LIMITS['beam_width']: return complete
                    elif len(next_states)<LIMITS['beam_width']:
                        next_states.append(changed.GetMol())
        if complete: break
        pending=next_states
        if not pending: break
    return complete

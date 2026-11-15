"""Finite graph interventions from explicit fixed identity and spatial adjacency.

No original scaffold is copied. Atom count, order, elements and coordinates are
preserved. Unknown generated-region bond orders and missing atoms are not guessed.
"""
from itertools import combinations
import numpy as np
from rdkit import Chem, rdBase
from .io import require

LIMITS = dict(max_candidates=4, max_sanitizations=128, max_fixed_external_removals=2,
              max_new_bonds=3, beam_width=4, max_bridge_pairs=24,
              max_bridge_distance_A=2.5, max_bridge_radius_ratio=1.6,
              min_bridge_radius_ratio=0.6)


def bond_map(mol):
    return {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): str(b.GetBondType()) for b in mol.GetBonds()}


def graph_key(mol):
    return (tuple((a.GetAtomicNum(), a.GetFormalCharge()) for a in mol.GetAtoms()), tuple(sorted(bond_map(mol).items())))


def graph_edits(original, changed):
    # Do not count Kekule/aromatic serialization differences as bond repairs.
    original = Chem.Mol(original)
    try:
        with rdBase.BlockLogs():
            Chem.SanitizeMol(original)
    except (ValueError, RuntimeError):
        pass
    a, b = bond_map(original), bond_map(changed)
    bonds = [dict(atoms=list(pair), before=a.get(pair), after=b.get(pair))
             for pair in sorted(a.keys() | b.keys()) if a.get(pair) != b.get(pair)]
    charges = [dict(atom=i, before=x.GetFormalCharge(), after=y.GetFormalCharge())
               for i, (x, y) in enumerate(zip(original.GetAtoms(), changed.GetAtoms()))
               if x.GetFormalCharge() != y.GetFormalCharge()]
    return dict(bonds=bonds, charges=charges)


def local_fixed(item):
    slots = {s: i for i, s in enumerate(item['slots'])}
    task = item['task']
    require(set(task['fixed_sampler_indices']) <= slots.keys(), 'Mapped fixed atoms are missing')
    return {slots[s]: r for s, r in zip(task['fixed_sampler_indices'], task['fixed_reference_indices'])}


def fixed_matches(mol, ref, mapping):
    for i, r in mapping.items():
        a, b = mol.GetAtomWithIdx(i), ref.GetAtomWithIdx(r)
        if (a.GetAtomicNum(), a.GetFormalCharge()) != (b.GetAtomicNum(), b.GetFormalCharge()):
            return False
    a, b = bond_map(mol), bond_map(ref)
    return all(a.get(tuple(sorted((i, j)))) == b.get(tuple(sorted((mapping[i], mapping[j]))))
               for i, j in combinations(mapping, 2))


class Budget:
    def __init__(self):
        self.attempts = 0
        self.rejected = 0

    def sanitize(self, mol):
        if self.attempts >= LIMITS['max_sanitizations']:
            return None
        self.attempts += 1
        candidate = Chem.Mol(mol)
        try:
            with rdBase.BlockLogs():
                candidate.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(candidate)
            return candidate
        except (ValueError, RuntimeError):
            self.rejected += 1
            return None


def restore_fixed_graph(mol, ref, mapping, budget):
    changed = Chem.RWMol(mol)
    for i, r in mapping.items():
        require(mol.GetAtomWithIdx(i).GetAtomicNum() == ref.GetAtomWithIdx(r).GetAtomicNum(),
                'Element mismatch cannot be repaired by graph editing')
        atom, target = changed.GetAtomWithIdx(i), ref.GetAtomWithIdx(r)
        atom.SetFormalCharge(target.GetFormalCharge())
        atom.SetIsAromatic(target.GetIsAromatic())
    for i, j in combinations(mapping, 2):
        old = changed.GetBondBetweenAtoms(i, j)
        target = ref.GetBondBetweenAtoms(mapping[i], mapping[j])
        if old is not None:
            changed.RemoveBond(i, j)
        if target is not None:
            changed.AddBond(i, j, target.GetBondType())
            changed.GetBondBetweenAtoms(i, j).SetIsAromatic(target.GetIsAromatic())
    first = budget.sanitize(changed)
    if first is not None and fixed_matches(first, ref, mapping):
        return [first]
    # Template restoration can reveal an incompatible extra external bond
    # (e.g. a fixed halogen was reattached elsewhere). Try minimal deletions.
    edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in changed.GetBonds()
             if (b.GetBeginAtomIdx() in mapping) != (b.GetEndAtomIdx() in mapping)]
    xyz = mol.GetConformer().GetPositions()
    edges.sort(key=lambda e: (-float(np.linalg.norm(xyz[e[0]]-xyz[e[1]])), tuple(sorted(e))))
    results = []
    for count in range(1, LIMITS['max_fixed_external_removals']+1):
        for removed in combinations(edges, count):
            if budget.attempts >= LIMITS['max_sanitizations']:
                return results
            trial = Chem.RWMol(changed)
            for a, b in removed:
                trial.RemoveBond(a, b)
            candidate = budget.sanitize(trial)
            if candidate is not None and fixed_matches(candidate, ref, mapping):
                results.append(candidate)
                if len(results) >= LIMITS['beam_width']:
                    return results
        if results:
            break
    return results


def bridge_pairs(mol, fixed):
    fragments = Chem.GetMolFrags(mol)
    groups = {i: group for group, atoms in enumerate(fragments) for i in atoms}
    xyz = mol.GetConformer().GetPositions()
    table = Chem.GetPeriodicTable()
    result = []
    for a, b in combinations(range(mol.GetNumAtoms()), 2):
        if groups[a] == groups[b] or (a in fixed and b in fixed):
            continue
        atom_a, atom_b = mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(b)
        # Covalent main-group graph repair; metal coordination is a separate method.
        if atom_a.GetSymbol() not in {'B','C','N','O','P','S','F','Cl','Br','I','Si'} or atom_b.GetSymbol() not in {'B','C','N','O','P','S','F','Cl','Br','I','Si'}:
            continue
        distance = float(np.linalg.norm(xyz[a]-xyz[b]))
        radius = table.GetRcovalent(atom_a.GetAtomicNum()) + table.GetRcovalent(atom_b.GetAtomicNum())
        ratio = distance / radius
        if (distance <= LIMITS['max_bridge_distance_A'] and
                LIMITS['min_bridge_radius_ratio'] <= ratio <= LIMITS['max_bridge_radius_ratio']):
            result.append((abs(ratio-1.), a, b, distance))
    return sorted(result)[:LIMITS['max_bridge_pairs']]


def reconnect(mol, fixed, budget):
    n = len(Chem.GetMolFrags(mol))
    if n <= 1:
        return []
    if n-1 > LIMITS['max_new_bonds']:
        return []
    beam = [(0., mol)]
    seen = {graph_key(mol)}
    completed = []
    for _ in range(n-1):
        next_states = []
        for score, base in beam:
            for penalty, a, b, _distance in bridge_pairs(base, fixed):
                if budget.attempts >= LIMITS['max_sanitizations']:
                    return completed
                changed = Chem.RWMol(base)
                changed.AddBond(a, b, Chem.BondType.SINGLE)
                candidate = budget.sanitize(changed)
                if candidate is None:
                    continue
                key = graph_key(candidate)
                if key in seen:
                    continue
                seen.add(key)
                next_states.append((score+penalty, candidate))
        next_states.sort(key=lambda v: (v[0], graph_key(v[1])))
        beam = next_states[:LIMITS['beam_width']]
        completed = [m for _, m in beam if len(Chem.GetMolFrags(m)) == 1]
        if completed or not beam:
            break
    return completed


def plan(item):
    mol, ref = item['mol'], item['ref']
    audit = dict(settings=LIMITS.copy(), attempted_sanitizations=0, rejected_graphs=0,
                 candidates=[], limitations=['No atom insertion/deletion or element replacement',
                    'No arbitrary generated-core bond-order reconstruction',
                    'No neural or trajectory intervention', 'Bridge proposals are hypotheses, not inferred true bonds'])
    if mol is None:
        audit['status'] = 'unreadable_molecule'
        return [], audit
    mapping = local_fixed(item)
    if any(mol.GetAtomWithIdx(i).GetAtomicNum() != ref.GetAtomWithIdx(r).GetAtomicNum() for i,r in mapping.items()):
        audit['status'] = 'fixed_element_mismatch'
        return [], audit
    budget = Budget()
    valid_original = budget.sanitize(mol)
    if valid_original is not None:
        mol = valid_original
    families = {'F': [], 'L': [], 'FL': [], 'E': [], 'EL': []}
    if not fixed_matches(mol, ref, mapping):
        families['F'] = restore_fixed_graph(mol, ref, mapping, budget)
    if valid_original is not None:
        families['L'] = reconnect(valid_original, set(mapping), budget)
    for candidate in families['F']:
        families['FL'].extend(reconnect(candidate, set(mapping), budget))
    if valid_original is None:
        from .chemistry import repair_invalid
        families['E'] = repair_invalid(mol,ref,mapping,budget)
        for candidate in families['E']:
            families['EL'].extend(reconnect(candidate,set(mapping),budget))
    # Reserve representation for each applicable mechanism; remaining slots use
    # edit count / graph order, never the final evaluator's outcome.
    selected, seen = [], {graph_key(mol)}
    queues = {k: sorted(v, key=lambda m: (len(graph_edits(mol,m)['bonds']), graph_key(m))) for k,v in families.items()}
    while len(selected) < LIMITS['max_candidates'] and any(queues.values()):
        for family in ('F', 'L', 'FL', 'E', 'EL'):
            while queues[family]:
                candidate = queues[family].pop(0)
                key = graph_key(candidate)
                if key in seen:
                    continue
                seen.add(key)
                name = f'{family}{1+sum(c["family"] == family for c in selected)}'
                require(candidate.GetNumAtoms() == mol.GetNumAtoms(), 'Topology changed atom count')
                require(np.array_equal(candidate.GetConformer().GetPositions(), mol.GetConformer().GetPositions()),
                        'Graph proposal changed coordinates')
                entry = dict(name=name, family=family, mol=candidate, edits=graph_edits(mol,candidate),
                             disconnected_components=len(Chem.GetMolFrags(candidate)))
                selected.append(entry)
                break
            if len(selected) >= LIMITS['max_candidates']:
                break
    audit.update(attempted_sanitizations=budget.attempts, rejected_graphs=budget.rejected,
                 initial_components=len(Chem.GetMolFrags(mol)),
                 initial_eligible_bridge_pairs=len(bridge_pairs(mol,set(mapping))),
                 prepared_by_family={k:len(v) for k,v in families.items()},
                 budget_exhausted=budget.attempts >= LIMITS['max_sanitizations'],
                 status='candidates_prepared' if selected else 'no_admissible_candidate_within_budget',
                 candidates=[{k:v for k,v in c.items() if k != 'mol'} for c in selected])
    return selected, audit

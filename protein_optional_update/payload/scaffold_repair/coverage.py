"""Literature-derived coverage, separated from the unchanged five-part endpoint.

Extra task/property metrics are not silently substituted for structural validity.
Unassessed properties are never reported as passes.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import QED, Crippen, Descriptors, rdFingerprintGenerator
from .chemistry import problems

SOURCES = {
 'DiffSBDD':'https://www.nature.com/articles/s43588-024-00737-x',
 'DiffHopp':'https://icml-compbio.github.io/2023/papers/WCBICML2023_paper69.pdf',
 'TurboHopp':'https://arxiv.org/abs/2410.20660',
 'PoseBusters':'https://pubs.rsc.org/en/content/articlehtml/2024/sc/d3sc04185a',
 'PoseCheck':'https://openreview.net/pdf/afaf5d2ccada64db4c2fcd8f0434ac552e1d563d.pdf',
 'FLOWR.ROOT':'https://www.nature.com/articles/s41467-026-74130-9',
 'MolGenBench':'https://www.nature.com/articles/s41467-026-78105-8',
 'RDKit':'https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html',
}

CATALOG = [
 dict(id='input_integrity', problems='Unreadable SDF, missing files/mapping, invalid coordinate frame',
      detection='Input validation; missing evidence stops, unreadable molecule retained',
      interventions=[], scope='Repair input/export upstream; molecular edits cannot reconstruct absent evidence',sources=['PoseBusters']),
 dict(id='chemical_graph', problems='Valence, bond order, aromaticity, fixed elements/charges',
      detection='RDKit chemistry problems with atom IDs; mapped fixed graph comparison',
      interventions=['F','E','FL','EL'],scope='Bounded template restoration or local invalid-bond hypotheses; no arbitrary generated-core reconstruction',sources=['PoseBusters','DiffSBDD']),
 dict(id='connectivity',problems='Disconnected components, missing or incompatible attachments',
      detection='Connected components; fixed-boundary signatures; spatial candidate pairs',
      interventions=['L','FL','EL'],scope='Nearby admissible single bonds only; no long missing-linker invention',sources=['DiffHopp','TurboHopp']),
 dict(id='intramolecular_geometry',problems='Bond length/angle, internal clash, aromatic/double-bond planarity, strain',
      detection='All configured PoseBusters intramolecular checks, local geometry details',
      interventions=['A','B','AB','C','G'],scope='Bounded force field relaxation or anchored conformer reconstruction',sources=['PoseBusters','PoseCheck']),
 dict(id='environment',problems='Protein, water, organic/inorganic cofactor clashes and volume overlap; wrong pocket placement',
      detection='All configured PoseBusters environment checks; residue-resolved distance pairs',
      interventions=['B','AB','C','G'],scope='No protein movement or deletion of real cofactors; frame mismatch requires corrected inputs',sources=['PoseBusters','PoseCheck']),
 dict(id='task_constraints',problems='Fixed-fragment drift, fixed stereochemistry, attachment signature, unchanged scaffold',
      detection='Mapping, displacement, local fully mapped tetrahedral geometry, boundary signatures, Murcko identity',
      interventions=['A','F','L','G'],scope='Only applicable structural interventions. Extra task flags remain separate; novel scaffold design is not a coordinate repair',sources=['DiffHopp','FLOWR.ROOT']),
 dict(id='chemical_properties',problems='Low drug-likeness, difficult synthesis, structural alerts',
      detection='QED, SA when available, logP, molecular weight, PAINS alerts',
      interventions=[],scope='Descriptive scores/alerts, not universal binary failures; no automatic efficacy/toxicity/synthesis claim',sources=['DiffSBDD','TurboHopp','MolGenBench']),
 dict(id='distribution',problems='Duplicates, low diversity, lack of novelty',
      detection='Batch uniqueness/diversity. Training-set novelty unassessed without comparison data',
      interventions=[],scope='Resampling/retraining needs generator access; deleting duplicates is not repairing molecules',sources=['DiffSBDD','TurboHopp']),
 dict(id='binding',problems='Poor binding/interaction retention, docking-workflow inconsistency, activity/selectivity',
      detection='Explicitly unassessed: no docking/interaction/experimental oracle configured',
      interventions=[],scope='Requires target-specific oracle and fixed protocol; structural pass does not establish activity',sources=['PoseCheck','MolGenBench']),
]


def geometry_details(mol, config):
    from posebusters.modules.distance_geometry import check_geometry
    from posebusters.modules.flatness import check_flatness
    details=[]
    for module in config['modules']:
        params=module.get('parameters',{})
        if module['function']=='distance_geometry':
            result=check_geometry(mol,**params)
            for kind,frame in result.get('details',{}).items():
                field,threshold={'bonds':('percent_error',params.get('threshold_bad_bond_length',.2)),
                    'angles':('bound_absolute_percent_error',params.get('threshold_bad_angle',.2)),
                    'clash':('bound_percent_error',params.get('threshold_clash',.2))}[kind]
                for row in frame.to_dict('records'):
                    value=float(row[field])
                    failed=value < -threshold if kind=='clash' else abs(value)>threshold
                    if failed:
                        atoms=row['angle'] if kind=='angles' else row['atom_pair']
                        details.append(dict(type=kind,atoms=list(map(int,atoms)),distance_A=float(row['distance']),
                            lower_bound_A=float(row['lower_bound']),upper_bound_A=float(row['upper_bound']),
                            observed_relative_error=value,threshold=threshold))
        elif module['function']=='flatness':
            d=check_flatness(mol,**params).get('details',{})
            for typ,atoms,distance,passed in zip(d.get('type',[]),d.get('planar_group',[]),d.get('max_distance',[]),d.get('flatness_passes',[])):
                if not passed:
                    details.append(dict(type=typ,atoms=list(map(int,atoms)),maximum_plane_distance_A=float(distance),
                        threshold_A=params['threshold_flatness']))
    return details


def fixed_stereo(item,mol):
    """Compare orientation only where center and its complete heavy neighborhood are mapped.

    CIP labels alone are not compared: priorities can change across a new scaffold.
    No claim is made for E/Z or partly mapped stereocenters.
    """
    from .topology import local_fixed
    mapping=local_fixed(item); inverse={r:i for i,r in mapping.items()}
    result=[]
    ref=item['ref']; x=mol.GetConformer().GetPositions(); y=ref.GetConformer().GetPositions()
    for i,r in mapping.items():
        atom=ref.GetAtomWithIdx(r)
        if atom.GetChiralTag()==Chem.ChiralType.CHI_UNSPECIFIED: continue
        neighbors=sorted(n.GetIdx() for n in atom.GetNeighbors())
        expected={inverse[n] for n in neighbors if n in inverse}
        actual={n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()}
        if not set(neighbors)<=inverse.keys() or actual!=expected or len(neighbors) not in (3,4):
            result.append(dict(atom=i,reference_atom=r,status='unassessed_incomplete_local_mapping'))
            continue
        a=x[[inverse[n] for n in neighbors]]; b=y[neighbors]
        v=float(np.linalg.det(a[:3]-(a[3] if len(a)==4 else x[i])))
        w=float(np.linalg.det(b[:3]-(b[3] if len(b)==4 else y[r])))
        status='degenerate_geometry' if abs(v)<1e-5 or abs(w)<1e-5 else ('same_orientation' if v*w>0 else 'inverted_orientation')
        result.append(dict(atom=i,reference_atom=r,status=status,generated_signed_volume=v,reference_signed_volume=w))
    return dict(centers=result,scope='Fully mapped tetrahedral neighborhoods only; E/Z and boundary stereo unassessed',
                automatically_repaired=False)


@lru_cache(maxsize=8192)
def properties(smiles):
    mol=Chem.MolFromSmiles(smiles)
    out=dict(qed=float(QED.qed(mol)),logp=float(Crippen.MolLogP(mol)),molecular_weight=float(Descriptors.MolWt(mol)))
    try:
        from rdkit.Contrib.SA_Score import sascorer
        out['sa_score']=float(sascorer.calculateScore(mol))
    except Exception as exc:
        out['sa_score']=None;out['sa_status']='unassessed: '+str(exc)
    try:
        from rdkit.Chem import FilterCatalog
        p=FilterCatalog.FilterCatalogParams();p.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        matches=FilterCatalog.FilterCatalog(p).GetMatches(mol)
        out['pains_alerts']=sorted({m.GetDescription() for m in matches})
    except Exception as exc:
        out['pains_alerts']=None;out['pains_status']='unassessed: '+str(exc)
    out['interpretation']='Descriptive only; PAINS alerts are not proof of toxicity, SA is not a synthesis route, QED is not activity'
    return out


def assess(item,mol,row,config):
    report=dict(structural_endpoint_pass=all(row.get(k) is True for k in
        ('rdkit_valid','fc_valid','root_fixed_substructure_match','fixed_chemistry_preserved','posebusters_valid')),
        chemistry_problems=problems(mol),geometry_locations=[],task_checks={},properties=None,
        binding=dict(status='not_assessed',reason='No docking, interaction or activity oracle configured'),
        inference_time_cause='unavailable_from_final_SDF',
        scope='Supplementary task/property coverage; not silently added to the five-part endpoint')
    report.update(evaluation_scope=row.get('evaluation_scope'), success_meaning=row.get('success_meaning'),
                  not_evaluated_checks=row.get('not_evaluated_checks', {}))
    if mol is not None and row['rdkit_valid']:
        report['geometry_locations']=geometry_details(mol,config)
        report['task_checks']=dict(scaffold_changed=row['murcko_scaffold_changed'],
            strict_attachment_signature_matches_reference=row.get('strict_boundary_reference_match'),
            attachment_details=row.get('boundary_details'),fixed_max_displacement_A=row.get('fixed_max_displacement_A'),
            fixed_stereochemistry=fixed_stereo(item,mol),
            interpretation='Attachment signatures and novel scaffold requirements depend on task; no universal threshold is invented')
        report['properties']=properties(Chem.MolToSmiles(mol,isomericSmiles=True))
    return report


def batch_metrics(smiles, total):
    valid=[s for s in smiles if s is not None]
    counts=Counter(valid)
    out=dict(total_records=total,rdkit_valid_records=len(valid),unique_molecules=len(counts),
        uniqueness_among_valid=len(counts)/len(valid) if valid else None,
        duplicate_counts={s:n for s,n in counts.items() if n>1},
        training_novelty=dict(status='not_assessed',reason='No training-set comparison database supplied'))
    if len(valid)>1:
        gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in valid]
        n=len(fps); total_pairs=n*(n-1)//2
        if total_pairs<=10000:
            pairs=list(combinations(range(n),2));estimated=False
        else:
            rng=np.random.default_rng(12037);selected=set()
            while len(selected)<10000:
                a,b=map(int,rng.integers(0,n,2))
                if a!=b: selected.add(tuple(sorted((a,b))))
            pairs=sorted(selected);estimated=True
        out['diversity_1_minus_morgan_tanimoto']=float(np.mean([1-DataStructs.TanimotoSimilarity(fps[a],fps[b]) for a,b in pairs]))
        out.update(diversity_pairs=len(pairs),diversity_estimated=estimated)
    else: out['diversity_1_minus_morgan_tanimoto']=None
    return out

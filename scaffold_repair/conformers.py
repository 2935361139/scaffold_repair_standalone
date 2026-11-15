"""Finite anchored conformer proposals; independent of the generator/model name."""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D
from .topology import graph_edits, fixed_matches, local_fixed, graph_key

SETTINGS=dict(max_candidates=2, seeds=[17011,29009], max_attempts=30,
              timeout_seconds=20, maximum_anchor_fit_rmsd_A=.2)


def propose(item, graph_candidates):
    audit=dict(settings=SETTINGS.copy(),attempts=[],status='no_eligible_connected_graph')
    mapping=local_fixed(item)
    eligible=[]
    for family,mol in [('baseline',item['mol'])]+[(c['name'],c['mol']) for c in graph_candidates]:
        if mol is None: continue
        mol=Chem.Mol(mol)
        try:
            with rdBase.BlockLogs(): Chem.SanitizeMol(mol)
        except (ValueError,RuntimeError): continue
        if len(Chem.GetMolFrags(mol))==1 and fixed_matches(mol,item['ref'],mapping):
            eligible.append((family,mol))
    if not eligible: return [],audit
    # Choose by predeclared priority, not final evaluator outcomes.
    source,base=eligible[0]
    original=item['xyz'][item['slots']].astype(float)
    fixed=sorted(mapping)
    if len(fixed)==len(original):
        audit['status']='no_variable_atoms'
        return [],audit
    proposals=[]
    for seed in SETTINGS['seeds'][:SETTINGS['max_candidates']]:
        record=dict(seed=seed,source_graph=source,status='not_embedded')
        audit['attempts'].append(record)
        try:
            heavy=Chem.Mol(base)
            heavy.RemoveAllConformers()
            hydrogenated=Chem.AddHs(heavy)
            params=AllChem.ETKDGv3()
            params.randomSeed=seed
            params.maxIterations=SETTINGS['max_attempts']
            params.timeout=SETTINGS['timeout_seconds']
            params.numThreads=1
            params.enforceChirality=True
            params.trackFailures=True
            params.SetCoordMap({i:Point3D(*map(float,original[i])) for i in fixed})
            with rdBase.BlockLogs(): cid=AllChem.EmbedMolecule(hydrogenated,params)
            record['embedding_failures']=list(params.GetFailureCounts())
            if cid<0:
                record['status']='embedding_failed_within_budget'
                continue
            candidate=Chem.Mol(base)
            conf=Chem.Conformer(len(original))
            for i in range(len(original)): conf.SetAtomPosition(i,hydrogenated.GetConformer().GetAtomPosition(i))
            candidate.RemoveAllConformers(); candidate.AddConformer(conf)
            reference=Chem.Mol(base)
            for i,x in enumerate(original): reference.GetConformer().SetAtomPosition(i,tuple(map(float,x)))
            rms=rdMolAlign.AlignMol(candidate,reference,atomMap=[(i,i) for i in fixed],reflect=False)
            record['anchor_fit_rmsd_A']=float(rms)
            if rms>SETTINGS['maximum_anchor_fit_rmsd_A']:
                record['status']='anchor_geometry_incompatible'
                continue
            xyz=candidate.GetConformer().GetPositions()
            xyz[fixed]=original[fixed]
            for i,x in enumerate(xyz): candidate.GetConformer().SetAtomPosition(i,tuple(map(float,x)))
            if graph_key(candidate)!=graph_key(base): raise ValueError('Embedding changed graph')
            name=('G' if source=='baseline' else source+'G')+str(len(proposals)+1)
            record['status']='candidate_prepared'; record['name']=name
            proposals.append(dict(name=name,family='G',mol=candidate,coordinates=xyz,
                edits=graph_edits(item['mol'],candidate),source_graph=source,
                disconnected_components=1,receipt=record))
        except Exception as exc:
            record.update(status='proposal_error',error=str(exc))
    audit['status']='candidates_prepared' if proposals else 'no_conformer_within_budget'
    return proposals,audit

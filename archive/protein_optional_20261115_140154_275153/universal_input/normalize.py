import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import numpy as np
from rdkit import Chem
from scaffold_repair.import_state import restore, block, bond_signature

BONDS = {'SINGLE': Chem.BondType.SINGLE, 'DOUBLE': Chem.BondType.DOUBLE,
         'TRIPLE': Chem.BondType.TRIPLE, 'AROMATIC': Chem.BondType.AROMATIC}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_json(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def write_json(p, data):
    Path(p).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def resolve(p, base):
    return (Path(base) / Path(p).expanduser()).resolve()


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_record(path, index=0):
    supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
    require(type(index) is int and 0 <= index < len(supplier), 'SDF record index out of range')
    return restore(supplier[index]), supplier.GetItemText(index)


def build_graph(a, decoder, bond_decoder, encoding):
    x, q, e, pos = [a[k] for k in ('atoms', 'charges', 'bonds', 'positions')]
    n = len(x)
    require(x.shape == (n,) and q.shape == (n,) and e.shape == (n,n) and pos.shape == (n,3), 'Bad graph array shapes')
    require(all(np.isfinite(v).all() for v in (x,q,e,pos)), 'Nonfinite graph arrays')
    require(all(np.array_equal(v, np.round(v)) for v in (x,q,e)), 'Graph labels must be integer categories')
    require(np.array_equal(e,e.T) and np.all(np.diag(e)==0), 'Bond matrix must be symmetric with zero diagonal')
    mol = Chem.RWMol()
    for atom, charge in zip(x,q):
        if encoding == 'index':
            require(0 <= int(atom) < len(decoder), 'Unknown atom category')
            at = Chem.Atom(decoder[int(atom)])
        else:
            require(encoding == 'atomic_number' and 1 <= int(atom) <= 118, 'Invalid atomic number/encoding')
            at = Chem.Atom(int(atom))
        at.SetFormalCharge(int(charge))
        mol.AddAtom(at)
    for i in range(n):
        for j in range(i+1,n):
            if int(e[i,j]) == 0:
                continue
            kind = bond_decoder.get(str(int(e[i,j])))
            require(kind in BONDS, 'Unknown bond category: ' + str(e[i,j]))
            mol.AddBond(i,j,BONDS[kind])
            if kind == 'AROMATIC':
                mol.GetAtomWithIdx(i).SetIsAromatic(True)
                mol.GetAtomWithIdx(j).SetIsAromatic(True)
    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i, xyz in enumerate(pos):
        conf.SetAtomPosition(i, tuple(map(float,xyz)))
    mol.AddConformer(conf)
    return mol.GetMol()


def heavy_view(mol):
    # Removing ordinary explicit H atoms is representation conversion. Preserve
    # their counts so an overvalent carbon cannot become valid by dropping H.
    for atom in mol.GetAtoms():
        atom.SetIntProp('_source_slot', atom.GetIdx())
        if atom.GetAtomicNum() == 1:
            require(atom.GetFormalCharge()==0 and atom.GetIsotope()==0 and atom.GetDegree()==1,
                    'Hydrogen charge/isotope/connectivity not representable by v2 heavy-atom input')
            require(atom.GetNeighbors()[0].GetAtomicNum()>1 and
                    atom.GetBonds()[0].GetBondType()==Chem.BondType.SINGLE,
                    'Nonstandard hydrogen bond cannot be silently removed')
    result = Chem.RemoveHs(mol, implicitOnly=False, updateExplicitCount=True, sanitize=False)
    require(all(a.GetAtomicNum()>1 for a in result.GetAtoms()), 'Hydrogens remain after conversion')
    mapping = {a.GetIntProp('_source_slot'):a.GetIdx() for a in result.GetAtoms()}
    return result, mapping


def config_text(manifest, experiment):
    return ('INPUT_MODE = "manifest"\nMANIFEST_PATH = '+repr(str(manifest))+
            '\nOUTPUT_DIR = '+repr(str(experiment))+
            '\nINTERVENTION_INPUT_DIR = None\nFAILURE_IDS = None\n')


def convert(cfg, base):
    base = Path(base)
    fmt = cfg['format']
    source = resolve(cfg['source'],base)
    out = resolve(cfg['output'],base)
    require(not out.exists(), 'Output exists; choose a new output directory: '+str(out))
    require(source.exists(), 'Source not found: '+str(source))
    if fmt == 'flowr_records':
        require(source.is_dir(), 'FLOWR source must be a complete generated directory')
        out.mkdir(parents=True)
        text = ('INPUT_MODE = "flowr_records"\nFLOWR_GENERATED_DIR = '+repr(str(source))+
                '\nPROTEIN_PDB = ""\nOUTPUT_DIR = '+repr(str(resolve(cfg['experiment_output'],base)))+
                '\nINTERVENTION_INPUT_DIR = None\nFAILURE_IDS = None\n')
        (out/'run_config.py').write_text(text)
        report = dict(total_samples=None,status='delegated_to_existing_flowr_adapter',needs=[],report=str(out/'conversion_report.json'))
        write_json(out/'conversion_report.json',report)
        return report
    out.mkdir(parents=True)
    write_json(out/'conversion_settings.json', cfg)
    receipts, specs, failures = [], [], []
    pending = set()
    source_copies = {}

    def snapshot(path):
        path = Path(path).resolve()
        require(path.is_file(), 'Input file missing: '+str(path))
        key = str(path)
        if key not in source_copies:
            target = out/'source'/('%03d_' % len(source_copies)+path.name)
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(path,target)
            source_copies[key] = {'copy':str(target.relative_to(out)), 'sha256':digest(path)}
        return path

    arrays = None
    if fmt == 'diffshape_capture':
        metadata = read_json(snapshot(source/'raw_manifest.json'))
        require(metadata.get('schema')=='diffshape_raw_capture_v1', 'Unknown capture schema')
        reference = snapshot(source/'template_original.sdf')
        require(digest(reference)==metadata['template_sha256'],'Template hash mismatch')
        ref,_ = read_record(reference)
        require(ref is not None,'Unreadable reference')
        changed = metadata['change_atom_idx_zero_based']
        require(len(set(changed))==len(changed) and all(type(i)is int and 0<=i<ref.GetNumAtoms() for i in changed),'Bad changed atom indices')
        fixed = [i for i in range(ref.GetNumAtoms()) if i not in changed and ref.GetAtomWithIdx(i).GetAtomicNum()>1]
        rows = [{'id':r['sample_id'],'fixed_generated_indices':fixed,'fixed_reference_indices':fixed} for r in metadata['samples']]
        require(metadata['total_generated']==len(rows)==metadata['requested_samples'],'Incomplete generation batch')
        arrays = np.load(snapshot(source/'raw_graphs.npz'), allow_pickle=False)
        profile = dict(atom_decoder=metadata['atom_decoder'],atom_encoding='index',
                       keys={'atoms':'{id}__atom_types','charges':'{id}__charges','bonds':'{id}__bond_types','positions':'{id}__positions'},
                       bond_decoder={'1':'SINGLE','2':'DOUBLE','3':'TRIPLE','4':'AROMATIC'})
        for row in rows:
            row['reference_sdf']=str(reference)
        # No rdkit_constructed, author success/failure labels, or diagnoses forwarded.
    elif fmt == 'npz_graphs':
        arrays = np.load(snapshot(source),allow_pickle=False)
        rows = read_json(snapshot(resolve(cfg['mapping_file'],base)))['samples']
        profile = cfg['graph_profile']
    elif fmt == 'sdf':
        rows = read_json(snapshot(resolve(cfg['mapping_file'],base)))['samples']
        supplier = Chem.SDMolSupplier(str(source),sanitize=False,removeHs=False)
        require(sorted(r['sdf_index'] for r in rows)==list(range(len(supplier))), 'Mapping must cover every SDF record exactly once; no silent filtering')
        snapshot(source)
        rows=[dict(r,sdf=str(source)) for r in rows]
    elif fmt == 'manifest':
        data = read_json(snapshot(source))
        require(data['schema']=='scaffold_repair_manifest_v1','Unknown manifest')
        rows=data['samples']
        rows=[dict(r,**{k:str(resolve(r[k],source.parent)) for k in ['sdf','reference_sdf','protein_pdb','coordinates_npy'] if r.get(k)}) for r in rows]
    else:
        raise ValueError('Unknown format: '+fmt)
    require(rows,'Empty batch')
    ids=[r['id'] for r in rows]
    require(len(set(ids))==len(ids) and all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,179}',s) for s in ids),'Invalid or duplicate sample IDs')
    for row in rows:
        sid=row['id']; folder=out/'samples'/sid;folder.mkdir(parents=True)
        try:
            require(not any(row.get(k) for k in ('trajectory_npz','trajectory_metadata','target_coordinates_npy')), 'Use the existing v2 manifest reader directly for trajectories or custom targets; this converter does not drop them silently')
            reference=snapshot(resolve(row.get('reference_sdf') or cfg['reference_sdf'],base))
            ref,_=read_record(reference,row.get('reference_index',0))
            require(ref is not None,'Unreadable reference');Chem.SanitizeMol(ref)
            require(ref.GetNumConformers()==1,'Reference needs one 3D conformer')
            fixed=list(row['fixed_generated_indices']);ref_fixed=list(row['fixed_reference_indices'])
            require(len(fixed)==len(ref_fixed)>0 and len(set(fixed))==len(fixed) and len(set(ref_fixed))==len(ref_fixed),'Invalid fixed mapping')
            require(all(type(i)is int and i>=0 for i in fixed+ref_fixed),'Indices must be nonnegative integers')
            if arrays is not None:
                a={k:np.asarray(arrays[v.format(id=sid,index=row.get('graph_index',ids.index(sid)))]) for k,v in profile['keys'].items()}
                np.savez_compressed(folder/'original_graph.npz',**a)
                mol=build_graph(a,profile.get('atom_decoder',[]),profile['bond_decoder'],profile.get('atom_encoding','index'))
            else:
                path=snapshot(resolve(row['sdf'],base));mol,raw=read_record(path,row.get('sdf_index',0))
                (folder/'original_record.sdf').write_text(raw)
                require(mol is not None,'Unreadable SDF record; retained, not dropped')
                if row.get('coordinates_npy'):
                    xyz=np.load(snapshot(resolve(row['coordinates_npy'],base)),allow_pickle=False)
                    old=mol.GetConformer().GetPositions()
                    require(xyz.shape==old.shape and np.isfinite(xyz).all() and np.max(abs(xyz-old))<=.002,'Coordinate array inconsistent with SDF')
                    for i,p in enumerate(xyz):mol.GetConformer().SetAtomPosition(i,tuple(map(float,p)))
            require(mol.GetNumConformers()==1,'Generated molecule needs one conformer')
            require(all(i<mol.GetNumAtoms() for i in fixed) and all(i<ref.GetNumAtoms() for i in ref_fixed),'Mapping outside atom range')
            require(all(mol.GetAtomWithIdx(i).GetAtomicNum()>1 and ref.GetAtomWithIdx(j).GetAtomicNum()>1 for i,j in zip(fixed,ref_fixed)),'Fixed mapping must refer to heavy atoms')
            require(cfg.get('units','angstrom') in ['angstrom','nanometer'],'Unsupported coordinate unit')
            scale=10. if cfg.get('units')=='nanometer' else 1.
            xyz=np.array(mol.GetConformer().GetPositions())*scale
            target=ref.GetConformer().GetPositions()[ref_fixed]
            require(np.isfinite(xyz).all() and np.isfinite(target).all(),'Nonfinite coordinates')
            frame=row.get('frame',cfg.get('frame','unverified'))
            rotation=np.eye(3);translation=np.zeros(3)
            if frame=='fixed_centroid_translation':
                require(all(mol.GetAtomWithIdx(i).GetAtomicNum()==ref.GetAtomWithIdx(j).GetAtomicNum() for i,j in zip(fixed,ref_fixed)),'Cannot register altered fixed elements')
                translation=np.mean(target-xyz[fixed],axis=0)
                residual=np.linalg.norm(xyz[fixed]+translation-target,axis=1)
                require(float(max(residual))<=cfg.get('maximum_registration_residual_A',.01),'Fixed atoms do not support a translation-only registration')
            elif frame=='explicit_rigid_transform':
                tf=row.get('transform',cfg.get('transform',{}))
                rotation=np.asarray(tf['rotation'],float);translation=np.asarray(tf['translation_A'],float)
                require(rotation.shape==(3,3) and translation.shape==(3,) and np.isfinite(rotation).all() and np.isfinite(translation).all(),'Invalid transform')
                require(np.allclose(rotation.T@rotation,np.eye(3),atol=1e-8) and abs(np.linalg.det(rotation)-1)<1e-8,'Transform must be a proper rigid rotation')
            elif frame not in ['reference','unverified']:
                raise ValueError('Unknown frame mode')
            if frame=='unverified':pending.add('Verify coordinate registration to reference/protein')
            xyz=xyz@rotation.T+translation
            for i,p in enumerate(xyz):mol.GetConformer().SetAtomPosition(i,tuple(map(float,p)))
            mol,hmap=heavy_view(mol);ref,rmap=heavy_view(ref)
            mapped=[hmap[i] for i in fixed];refmapped=[rmap[i] for i in ref_fixed]
            # Update caches without sanitizing or repairing; preserve valence/H semantics.
            mol.UpdatePropertyCache(strict=False)
            ref.UpdatePropertyCache(strict=False)
            (folder/'molecule.sdf').write_text(block(mol))
            recovered,_=read_record(folder/'molecule.sdf')
            require(recovered is not None and recovered.GetNumAtoms()==mol.GetNumAtoms(),'SDF cannot represent this graph')
            recovered.UpdatePropertyCache(strict=False)
            for a,b in zip(mol.GetAtoms(),recovered.GetAtoms()):
                require((a.GetAtomicNum(),a.GetFormalCharge(),a.GetTotalNumHs())==(b.GetAtomicNum(),b.GetFormalCharge(),b.GetTotalNumHs()),'SDF atom/HCOUNT roundtrip mismatch')
            require(bond_signature(mol)==bond_signature(recovered),'SDF bond roundtrip mismatch')
            (folder/'reference.sdf').write_text(block(ref))
            ref_recovered,_=read_record(folder/'reference.sdf')
            require(ref_recovered is not None and bond_signature(ref)==bond_signature(ref_recovered),'Reference SDF roundtrip mismatch')
            Chem.SanitizeMol(ref_recovered)
            np.save(folder/'coordinates.npy',mol.GetConformer().GetPositions())
            protein=row.get('protein_pdb') or cfg.get('protein_pdb')
            if protein:
                pp=snapshot(resolve(protein,base));shutil.copy2(pp,folder/'protein.pdb')
                require(any(l.startswith('ATOM  ') for l in pp.read_text().splitlines()),'PDB has no protein ATOM records')
                if not cfg.get('protein_frame_confirmed',False):pending.add('Confirm protein and reference use the same frame')
            else:pending.add('Provide matching protein PDB; current v2 requires it')
            relative=lambda name:str((folder/name).relative_to(out))
            spec=dict(id=sid,input_id=row.get('input_id',cfg.get('input_id','external_task')),model=cfg.get('model','external'),
                      sdf=relative('molecule.sdf'),sdf_index=0,reference_sdf=relative('reference.sdf'),
                      protein_pdb=relative('protein.pdb'),coordinates_npy=relative('coordinates.npy'),
                      fixed_generated_indices=mapped,fixed_reference_indices=refmapped)
            specs.append(spec)
            receipt=dict(id=sid,status='converted',heavy_slot_map=hmap,reference_heavy_slot_map=rmap,
                         frame_operation=frame,rotation=rotation.tolist(),translation_A=translation.tolist(),unit_scale=scale,
                         registration_residual_max_A=float(max(np.linalg.norm(xyz[fixed]-target,axis=1))),
                         interpretation='Coordinate registration is preprocessing; no atoms snapped to reference; no quality label used.')
            write_json(folder/'conversion_receipt.json',receipt);receipts.append(receipt)
        except Exception as exc:
            error={'id':sid,'status':'not_representable_or_input_error','reason':str(exc)}
            failures.append(error);write_json(folder/'conversion_error.json',error)
    if arrays is not None:arrays.close()
    if failures:pending.add('Resolve conversion errors; no partial batch sent to v2')
    ready=not pending
    name='manifest.json' if ready else 'manifest.pending.json'
    write_json(out/name,{'schema':'scaffold_repair_manifest_v1','samples':specs})
    if ready:
        (out/'run_config.py').write_text(config_text(out/name,resolve(cfg['experiment_output'],base)))
    report=dict(total_samples=len(rows),converted_samples=len(specs),conversion_errors=failures,
                status='ready' if ready else 'needs_context_or_representation',needs=sorted(pending),
                source_files=source_copies,report=str(out/'conversion_report.json'),
                no_author_labels_used=True,no_molecular_repairs=True,
                interpretation='Conversion success is NOT molecular validity. Retain all inputs and report conversion failures separately.')
    write_json(out/'conversion_report.json',report)
    return report

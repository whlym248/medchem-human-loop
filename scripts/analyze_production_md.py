"""Analyze one declared MD replicate, respecting v2 committed segment prefixes.

Dependencies: OpenMM, NumPy, MDTraj, Matplotlib, RDKit. No GPU is needed.
Does not edit input systems, trajectories, manifests or production reviews.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import sys

import numpy as np
import mdtraj as md
import openmm as mm
from openmm import app, unit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import ChemicalFeatures

PROTEIN=set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP ASH GLH LYN CYX CYM'.split())


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode('utf-8')).hexdigest()


def safe_path(root,relative):
    candidate=(root/relative).resolve()
    if not candidate.is_relative_to(root.resolve()):
        raise ValueError('Segment path escapes run directory')
    return candidate


def csv_prefix(path,n):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        fields=next(csv.reader([stream.readline().lstrip('#')]))
        rows=list(itertools.islice(csv.DictReader(stream,fieldnames=fields),n))
    if len(rows)!=n:
        raise ValueError(f'Not enough committed CSV rows: {path}')
    steps=np.array([int(float(r['Step'])) for r in rows],dtype=np.int64)
    times=np.array([float(r['Time (ps)']) for r in rows],dtype=float)
    if not np.isfinite(times).all() or len(steps) and np.any(np.diff(steps)<=0):
        raise ValueError(f'Nonfinite time or non-monotonic CSV steps: {path}')
    return steps,times


def load_segments(root,run,source,args):
    inputs={name:sha(source/name) for name in ('system.xml','complex.cif','positions_nm.npy','metadata.json')}
    manifest_path=run/'segment_manifest.json'
    records=[]; snapshot={}; warnings=[]
    if manifest_path.is_file():
        raw=manifest_path.read_bytes()
        manifest=json.loads(raw)
        identity=read_json(run/'run_identity.json')
        if identity.get('input_sha256')!=inputs:
            raise ValueError('Input hashes differ from v2 run identity')
        if identity.get('case_id')!=args.case:
            raise ValueError('v2 case mismatch')
        if manifest.get('run_identity_sha256')!=digest(identity):
            raise ValueError('Segment manifest identity mismatch')
        summary=read_json(run/'run_summary.json') if (run/'run_summary.json').is_file() else {}
        if not args.allow_incomplete and summary.get('status')!='completed':
            raise ValueError('Run not completed: use --allow-incomplete for an explicitly provisional snapshot')
        if summary.get('status')=='completed' and (summary.get('run_identity_sha256')!=digest(identity) or summary.get('completed_absolute_step')!=identity['target_absolute_step']):
            raise ValueError('Completed summary does not match run identity/target')
        dt=float(identity['settings']['timestep_ps']); eq=int(identity['equilibration_steps'])
        target=int(identity['target_absolute_step'])
        mode=identity['mode']; seed=identity['seed']
        for seg in manifest['segments']:
            if seg.get('run_identity_sha256')!=digest(identity):
                raise ValueError('Segment entry identity mismatch')
            n=int(seg['include_frame_count'])
            if n<0 or int(seg.get('include_first_frame_index_zero_based',0))!=0:
                raise ValueError('Unsupported segment prefix specification')
            if not n:continue
            if n>int(seg.get('paired_frames_committed',n)):
                raise ValueError('Included count exceeds paired committed frames')
            csvfile=safe_path(run,seg['state_csv']); dcd=safe_path(run,seg['trajectory'])
            steps,times=csv_prefix(csvfile,n)
            interval=int(seg['report_interval_steps'])
            first=(int(seg['include_steps_gt'])//interval+1)*interval
            if not np.array_equal(steps,first+np.arange(n)*interval):
                raise ValueError('CSV prefix does not match manifest absolute-step sequence')
            if steps[-1]>int(seg['include_steps_le']) or steps[0]<=eq or steps[-1]>target:
                raise ValueError('CSV prefix outside committed/production interval')
            if not np.allclose(times,steps*dt,rtol=0,atol=1e-3):
                raise ValueError('CSV time differs from absolute step * timestep')
            records.append(dict(attempt=seg['attempt'],csv=csvfile,dcd=dcd,steps=steps,times=times,
                                count=n,interval=interval,ignored_superseded_tail=seg.get('superseded_tail_after_step')))
        snapshot={'segment_manifest_sha256':hashlib.sha256(raw).hexdigest(),'manifest':manifest,
                  'run_identity':identity,'run_identity_sha256':sha(run/'run_identity.json'),
                  'run_summary':summary,'legacy':False}
        if summary.get('status')!='completed':warnings.append('Incomplete live snapshot; later committed frames can replace existing tails. Reanalyze final manifest.')
    else:
        if not args.legacy_engineering:
            raise ValueError('No segment manifest: legacy output requires --legacy-engineering')
        status=read_json(run/'run_status.json')
        if status.get('case_id')!=args.case or status.get('mode')!='pilot' or status.get('status')!='completed':
            raise ValueError('Legacy fallback is restricted to a completed matching engineering pilot')
        if status.get('system_sha256')!=inputs['system.xml']:
            raise ValueError('Legacy system hash mismatch')
        dt=0.002; eq=round(float(status['equilibration_ps'])/dt)
        with (run/'state.csv').open(encoding='utf-8-sig') as f:n=sum(1 for line in f)-1
        steps,times=csv_prefix(run/'state.csv',n)
        if n<2:raise ValueError('At least two legacy engineering frames required')
        intervals=np.diff(steps)
        if not np.all(intervals==intervals[0]) or not np.allclose(times,steps*dt,atol=1e-3):
            raise ValueError('Unexpected legacy engineering sampling grid')
        records=[dict(attempt='legacy_engineering',csv=run/'state.csv',dcd=run/'trajectory.dcd',steps=steps,times=times,
                      count=n,interval=int(intervals[0]),ignored_superseded_tail=None)]
        target=int(steps[-1]);mode='pilot';seed=status['seed']
        snapshot={'legacy':True,'legacy_status':status,'legacy_status_sha256':sha(run/'run_status.json')}
        warnings.append('Legacy engineering output used for software validation only; not production evidence.')
    if not records:raise ValueError('No committed sampled frames')
    allsteps=np.concatenate([s['steps'] for s in records])
    if np.any(np.diff(allsteps)<=0):
        raise ValueError('Duplicate or backward steps across segments; refusing to concatenate')
    previous=eq; covered=0; gaps=[]
    for rec in records:
        for step in rec['steps']:
            start=max(eq,int(step)-rec['interval'])
            if start>previous:gaps.append({'start_step':previous,'end_step':start})
            if start<previous:raise ValueError('Overlapping retained sampling intervals')
            covered+=int(step)-start; previous=int(step)
    if gaps:warnings.append('Sampling gaps exist: occupancy refers to observed retained frames, not unobserved intervals.')
    details={'input_sha256':inputs,'mode':mode,'seed':seed,'timestep_ps':dt,'equilibration_steps':eq,
             'excluded_equilibration_ns':eq*dt/1000,'declared_sampling_ns':(target-eq)*dt/1000,
             'target_absolute_step':target,'retained_frames':len(allsteps),'covered_sampling_ns':covered*dt/1000,
             'last_saved_absolute_time_ns':int(allsteps[-1])*dt/1000,
             'first_saved_production_time_ns':(int(allsteps[0])-eq)*dt/1000,
             'last_saved_production_time_ns':(int(allsteps[-1])-eq)*dt/1000,
             'saved_frame_span_ns':(int(allsteps[-1])-int(allsteps[0]))*dt/1000,
             'sampled_interval_gaps':gaps,'warnings':warnings,'snapshot':snapshot}
    return records,details


def covalent_bonds(system):
    bonds=set()
    for i in range(system.getNumConstraints()):
        a,b,_=system.getConstraintParameters(i);bonds.add(tuple(sorted((int(a),int(b)))))
    for force in system.getForces():
        if isinstance(force,mm.HarmonicBondForce):
            for i in range(force.getNumBonds()):
                a,b,_,k=force.getBondParameters(i)
                if k.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2)>0:
                    bonds.add(tuple(sorted((int(a),int(b)))))
    return sorted(bonds)


def mic(delta,box):
    fractional=np.asarray(delta)@np.linalg.inv(box)
    centered=fractional-np.rint(fractional)
    best=centered@box
    if np.allclose(box@box.T,np.diag(np.diag(box@box.T)),atol=1e-8):return best
    norm=np.sum(best*best,axis=-1)
    for shift in itertools.product((-1,0,1),repeat=3):
        trial=(centered+shift)@box; other=np.sum(trial*trial,axis=-1)
        improve=other<norm;best=np.where(improve[...,None],trial,best);norm=np.minimum(norm,other)
    return best


def make_components(n,bonds):
    graph=[set() for _ in range(n)]
    for a,b in bonds:graph[a].add(b);graph[b].add(a)
    unseen=set(range(n));trees=[]
    while unseen:
        first=min(unseen);unseen.remove(first);queue=[first];edges=[]
        for a in queue:
            for b in sorted(graph[a]):
                if b in unseen:unseen.remove(b);queue.append(b);edges.append((a,b))
        trees.append((np.array(queue,dtype=int),edges))
    return graph,trees


def whole_coordinates(xyz,box,trees,protein,ligand,pocket,reference):
    result=xyz.copy()
    for ids,edges in trees:
        if edges:
            pairs=np.asarray(edges,dtype=int)
            offsets=mic(xyz[pairs[:,1]]-xyz[pairs[:,0]],box)
            for (a,b),delta in zip(edges,offsets):result[b]=result[a]+delta
    pset=set(protein);lset=set(ligand)
    pcomponents=[ids for ids,_ in trees if any(i in pset for i in ids)]
    anchor=max(pcomponents,key=len); anchor_center=result[anchor].mean(axis=0)
    for ids in pcomponents:
        if np.array_equal(ids,anchor):continue
        center=result[ids].mean(axis=0)
        desired=anchor_center+reference[ids].mean(axis=0)-reference[anchor].mean(axis=0)
        result[ids]+=desired+mic(center-desired,box)-center
    pocket_center=result[pocket].mean(axis=0)
    for ids,_ in trees:
        if any(i in lset for i in ids):
            center=result[ids].mean(axis=0)
            result[ids]+=pocket_center+mic(center-pocket_center,box)-center
    return result


def label(atom):
    return f'{atom.residue.chain.id}:{atom.residue.name}{atom.residue.id}:{atom.name}[{atom.index}]'


def residue_label(residue):
    return f'{residue.chain.id}:{residue.name}{residue.id}'


def ligand_roles(source_json,all_atoms,ligand_ids):
    data=read_json(source_json)
    lookup={all_atoms[i].name:i for i in ligand_ids}
    if len(lookup)!=len(ligand_ids) or len(data['atoms'])!=len(ligand_ids):raise ValueError('Ambiguous ligand atom names')
    molecule=Chem.RWMol();mapping=[]
    for entry in data['atoms']:
        atom=Chem.Atom(int(entry['atomic_number']));atom.SetFormalCharge(int(entry['formal_charge']))
        atom.SetIsAromatic(bool(entry['is_aromatic']));molecule.AddAtom(atom)
        g=lookup[entry['name']]
        if all_atoms[g].element.atomic_number!=entry['atomic_number']:raise ValueError('Ligand element mapping mismatch')
        mapping.append(g)
    orders={1:Chem.BondType.SINGLE,2:Chem.BondType.DOUBLE,3:Chem.BondType.TRIPLE}
    for bond in data['bonds']:
        kind=Chem.BondType.AROMATIC if bond['is_aromatic'] else orders[bond['bond_order']]
        molecule.AddBond(int(bond['atom1']),int(bond['atom2']),kind)
    molecule=molecule.GetMol();Chem.SanitizeMol(molecule)
    factory=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
    donors=set();acceptors=set()
    for feature in factory.GetFeaturesForMol(molecule):
        if feature.GetFamily() in ('Donor','Acceptor'):
            target=donors if feature.GetFamily()=='Donor' else acceptors
            target.update(mapping[i] for i in feature.GetAtomIds() if molecule.GetAtomWithIdx(i).GetAtomicNum() in (7,8))
    return donors,acceptors


def hydrogen_bonds(all_atoms,selected,local,graph,ligand_global,ligand_donors,ligand_acceptors):
    ligand=set(ligand_global);protein_donors=set();protein_acceptors=set()
    def hs(g):return [selected[h] for h in graph[local[g]] if all_atoms[selected[h]].element.atomic_number==1]
    for g in selected:
        a=all_atoms[g]
        if g in ligand or a.residue.name not in PROTEIN:continue
        if a.element.atomic_number in (7,8) and hs(g):protein_donors.add(g)
        if a.element.atomic_number==8:
            acid=a.residue.name in ('ASP','ASH','GLU','GLH') and a.name not in ('O','OXT')
            if not ((acid or a.name=='OXT') and hs(g)):protein_acceptors.add(g)
        if a.residue.name in ('HIS','HID','HIE','HIP') and a.name in ('ND1','NE2') and not hs(g):
            protein_acceptors.add(g)
    groups=[];triplets=[];starts=[];triplet_group=[]
    for donors,acceptors,side in ((protein_donors,ligand_acceptors,'protein_to_ligand'),(ligand_donors,protein_acceptors,'ligand_to_protein')):
        for d in sorted(donors):
            hydrogens=hs(d)
            if not hydrogens:continue
            for a in sorted(acceptors):
                starts.append(len(triplets));gi=len(groups)
                pa=all_atoms[d if d not in ligand else a]
                groups.append({'donor_global':d,'acceptor_global':a,'donor':label(all_atoms[d]),'acceptor':label(all_atoms[a]),
                               'direction':side,'protein_residue':residue_label(pa.residue),'protein_resid':pa.residue.id})
                for h in sorted(hydrogens):triplets.append((local[d],local[h],local[a]));triplet_group.append(gi)
    return groups,np.asarray(triplets,dtype=int).reshape(-1,3),np.asarray(starts,dtype=int),np.asarray(triplet_group,dtype=int)


def hbond_geometry(chunk,group_pairs,triplets,group_starts,triplet_groups,args):
    """Return per-frame direct N/O geometry, including a valid zero-group case."""
    if not len(group_pairs):
        # Avoid empty MDTraj geometry calls and reduceat. No candidate pairs
        # means no observations; downstream occupancy tables stay empty.
        shape=(len(chunk),0)
        return (np.zeros(shape),np.zeros(shape),np.zeros(shape),
                np.zeros(shape,dtype=bool),np.zeros(shape,dtype=bool),np.zeros(shape))
    da=md.compute_distances(chunk,group_pairs,periodic=True)*10
    ha=md.compute_distances(chunk,triplets[:,[1,2]],periodic=True)*10
    angles=np.degrees(md.compute_angles(chunk,triplets,periodic=True))
    qualified=(da[:,triplet_groups]<=args.hbond_DA_A)&(ha<=args.hbond_HA_A)&(angles>=args.hbond_angle_deg)
    group_qualified=np.maximum.reduceat(qualified,group_starts,axis=1)
    group_angles=np.maximum.reduceat(angles,group_starts,axis=1)
    return da,ha,angles,qualified,group_qualified,group_angles


def write_table(path,rows,fields=None):
    if fields is None:fields=list(rows[0]) if rows else []
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def analyze(args):
    root=args.bundle_root.resolve();run=args.run_dir.resolve();out=args.out_dir.resolve()
    if '/' in args.case or '\\' in args.case or args.case in ('.','..'):raise ValueError('Invalid case')
    if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be new/empty to preserve prior analyses')
    out.mkdir(parents=True,exist_ok=True)
    source=root/'systems'/args.case;metadata=read_json(source/'metadata.json')
    records,details=load_segments(root,run,source,args)
    print(f'Committed frames to analyze: {details["retained_frames"]}',flush=True)
    pdb=app.PDBxFile(str(source/'complex.cif'));all_atoms=list(pdb.topology.atoms())
    system=mm.XmlSerializer.deserialize((source/'system.xml').read_text(encoding='utf-8'))
    if system.getNumParticles()!=len(all_atoms):raise ValueError('Atom count differs between system and topology')
    ligand_global=list(metadata['ligand_atom_indices']);ligand_set=set(ligand_global)
    selected=np.array([a.index for a in all_atoms if a.residue.name in PROTEIN or a.index in ligand_set],dtype=int)
    local={int(g):i for i,g in enumerate(selected)}
    protein=np.array([local[g] for g in selected if g not in ligand_set],dtype=int)
    ligand=np.array([local[g] for g in ligand_global],dtype=int)
    pheavy=np.array([i for i in protein if all_atoms[selected[i]].element.atomic_number>1],dtype=int)
    lheavy=np.array([i for i in ligand if all_atoms[selected[i]].element.atomic_number>1],dtype=int)
    backbone=np.array([i for i in protein if all_atoms[selected[i]].name in ('N','CA','C')],dtype=int)
    if len(backbone)<3 or not len(lheavy):raise ValueError('Insufficient backbone or ligand atoms')
    all_bonds=covalent_bonds(system)
    local_bonds=[(local[a],local[b]) for a,b in all_bonds if a in local and b in local]
    graph,trees=make_components(len(selected),local_bonds)
    mdtop=md.Topology.from_openmm(pdb.topology)
    existing={tuple(sorted((a.index,b.index))) for a,b in mdtop.bonds}
    for a,b in all_bonds:
        if (a,b) not in existing:mdtop.add_bond(mdtop.atom(a),mdtop.atom(b))
    reference=np.load(source/'positions_nm.npy',allow_pickle=False)[selected]
    refbox=np.asarray(pdb.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer),dtype=float)
    refdist=np.linalg.norm(mic(reference[pheavy,None,:]-reference[lheavy][None,:,:],refbox),axis=-1)
    pocket=pheavy[np.min(refdist,axis=1)<0.8]
    if len(pocket)<3:raise ValueError('No initial protein pocket near ligand')
    reference=whole_coordinates(reference,refbox,trees,protein,ligand,pocket,reference)
    refcenter=reference[backbone].mean(axis=0);y=reference[backbone]-refcenter
    residue_keys=[];residue_starts=[]
    for position,i in enumerate(pheavy):
        key=residue_label(all_atoms[selected[i]].residue)
        if not residue_keys or key!=residue_keys[-1]:residue_keys.append(key);residue_starts.append(position)
    if len(set(residue_keys))!=len(residue_keys):raise ValueError('Protein residue labels not unique')
    contact_pairs=np.array([(p,l) for p in pheavy for l in lheavy],dtype=int)
    initial_min=np.minimum.reduceat(np.min(refdist,axis=1),residue_starts)
    keyres=set(getattr(args,'highlight_residue',[]) or [])
    chosen=[i for i,name in enumerate(residue_keys) if initial_min[i]<=.5 or all_atoms[selected[pheavy[residue_starts[i]]]].residue.id in keyres]
    donors,acceptors=ligand_roles(root/'inputs'/'charged_ligands'/f'{args.case}.json',all_atoms,ligand_global)
    groups,triplets,group_starts,triplet_groups=hydrogen_bonds(all_atoms,selected,local,graph,ligand_global,donors,acceptors)
    group_pairs=np.array([(local[g['donor_global']],local[g['acceptor_global']]) for g in groups],dtype=int).reshape(-1,2)
    nh=len(groups);nres=len(residue_keys)
    ccount=np.zeros(nres,dtype=np.int64);csum=np.zeros(nres);cmin=np.full(nres,np.inf)
    hcount=np.zeros(nh,dtype=np.int64);hnear=np.zeros(nh,dtype=np.int64);hsum=np.zeros(nh);hmin=np.full(nh,np.inf);hangle=np.zeros(nh)
    rmsd_rows=[];contact_rows=[];done=0;input_records=[]
    hbfile=(out/'hydrogen_bond_geometry_near.csv').open('w',encoding='utf-8-sig',newline='')
    hbwriter=csv.writer(hbfile);hbwriter.writerow(['production_time_ns','absolute_step','segment','donor','hydrogen','acceptor','protein_residue','D_A_A','H_A_A','D_H_A_degrees','qualified'])
    try:
        for rec in records:
            dcd_size=rec['dcd'].stat().st_size
            used=0
            print(f'Analyzing {rec["attempt"]}: first {rec["count"]} paired frames',flush=True)
            for chunk in md.iterload(str(rec['dcd']),top=mdtop,atom_indices=selected,chunk=args.chunk_frames):
                take=min(len(chunk),rec['count']-used)
                if take<=0:break
                chunk=chunk[:take]
                if chunk.unitcell_vectors is None or not np.isfinite(chunk.xyz).all():raise ValueError('Missing periodic box or nonfinite DCD coordinates')
                steps=rec['steps'][used:used+take]
                times=(steps-details['equilibration_steps'])*details['timestep_ps']/1000
                distances=md.compute_distances(chunk,contact_pairs,periodic=True).reshape(take,len(pheavy),len(lheavy))
                residue_distance=np.minimum.reduceat(np.min(distances,axis=2),residue_starts,axis=1)
                ccount+=(residue_distance<=args.contact_A/10).sum(axis=0);csum+=residue_distance.sum(axis=0)*10;cmin=np.minimum(cmin,residue_distance.min(axis=0)*10)
                da,ha,angles,qualified,group_qualified,group_angles=hbond_geometry(chunk,group_pairs,triplets,group_starts,triplet_groups,args)
                hcount+=group_qualified.sum(axis=0);hnear+=(da<=5.0).sum(axis=0);hsum+=da.sum(axis=0);hmin=np.minimum(hmin,da.min(axis=0));hangle=np.maximum(hangle,group_angles.max(axis=0))
                for frame in range(take):
                    box=np.asarray(chunk.unitcell_vectors[frame],dtype=float)
                    whole=whole_coordinates(chunk.xyz[frame].astype(float),box,trees,protein,ligand,pocket,reference)
                    center=whole[backbone].mean(axis=0);x=whole[backbone]-center
                    u,s,vt=np.linalg.svd(x.T@y);fix=np.eye(3);fix[2,2]=np.linalg.det(u@vt)
                    aligned=(whole-center)@(u@fix@vt)+refcenter
                    rms=lambda ids:float(np.sqrt(np.mean(np.sum((aligned[ids]-reference[ids])**2,axis=1)))*10)
                    rmsd_rows.append({'absolute_step':int(steps[frame]),'production_time_ns':float(times[frame]),'segment':rec['attempt'],
                                      'protein_backbone_RMSD_A':rms(backbone),'ligand_heavy_RMSD_A':rms(lheavy),
                                      'protein_ligand_min_heavy_A':float(residue_distance[frame].min()*10)})
                    row={'absolute_step':int(steps[frame]),'production_time_ns':float(times[frame]),'segment':rec['attempt']}
                    row.update({residue_keys[r]+'_min_A':float(residue_distance[frame,r]*10) for r in chosen});contact_rows.append(row)
                    for t in np.flatnonzero(da[frame,triplet_groups]<=5.0):
                        g=int(triplet_groups[t]);d,h,a=triplets[t]
                        hbwriter.writerow([float(times[frame]),int(steps[frame]),rec['attempt'],label(all_atoms[selected[d]]),label(all_atoms[selected[h]]),label(all_atoms[selected[a]]),groups[g]['protein_residue'],float(da[frame,g]),float(ha[frame,t]),float(angles[frame,t]),int(qualified[frame,t])])
                used+=take;done+=take
                if used>=rec['count']:break
            if used!=rec['count']:raise ValueError('DCD has fewer frames than committed CSV prefix')
            input_records.append({'attempt':rec['attempt'],'dcd':str(rec['dcd']),'csv':str(rec['csv']),'retained_frames':used,
                                  'dcd_size_at_start':dcd_size,'dcd_size_at_end':rec['dcd'].stat().st_size,
                                  'retained_csv_prefix_sha256':hashlib.sha256((repr(rec['steps'].tolist())+'|'+repr(rec['times'].tolist())).encode()).hexdigest(),
                                  'ignored_superseded_tail_after_step':rec['ignored_superseded_tail']})
    finally:hbfile.close()
    if done!=details['retained_frames']:raise ValueError('Analyzed frame count mismatch')
    contacts=[{'residue':name,'initial_min_heavy_A':float(initial_min[i]*10),'observed_min_heavy_A':float(cmin[i]),
               'mean_min_heavy_A':float(csum[i]/done),'contact_frames':int(ccount[i]),'observed_frames':done,'contact_occupancy':float(ccount[i]/done)} for i,name in enumerate(residue_keys)]
    hbonds=[]
    for i,g in enumerate(groups):
        if hnear[i] or g['protein_resid'] in keyres:
            hbonds.append({**g,'qualified_frames':int(hcount[i]),'observed_frames':done,'occupancy':float(hcount[i]/done),
                           'frames_DA_le_5A':int(hnear[i]),'mean_DA_A':float(hsum[i]/done),'min_DA_A':float(hmin[i]),
                           'maximum_DHA_angle_any_donor_H_degrees':float(hangle[i])})
    hbonds.sort(key=lambda r:(-r['occupancy'],r['min_DA_A']))
    write_table(out/'rmsd_timeseries.csv',rmsd_rows)
    write_table(out/'key_contacts_timeseries.csv',contact_rows)
    write_table(out/'residue_contact_occupancy.csv',contacts)
    write_table(out/'hydrogen_bond_occupancy.csv',hbonds,fields=list(hbonds[0]) if hbonds else ['donor','acceptor','occupancy'])
    scope='ENGINEERING TEST' if details['mode']=='pilot' else ('PROVISIONAL MD SNAPSHOT' if args.allow_incomplete else 'MD REPLICATE')
    fig,axes=plt.subplots(2,1,figsize=(9,6),sharex=True,layout='constrained')
    for rec in records:
        rows=[r for r in rmsd_rows if r['segment']==rec['attempt']]
        axes[0].plot([r['production_time_ns'] for r in rows],[r['protein_backbone_RMSD_A'] for r in rows],color='#2667a4',linewidth=1.1)
        axes[1].plot([r['production_time_ns'] for r in rows],[r['ligand_heavy_RMSD_A'] for r in rows],color='#ce6a24',linewidth=1.1)
    axes[0].set_ylabel('Protein backbone RMSD (A)');axes[1].set_ylabel('Ligand heavy RMSD (A)');axes[1].set_xlabel('Unrestrained sampling time (ns)')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle(f'{scope} | {args.case} | seed {details["seed"]}\nPBC reconstructed; fit to input protein N/CA/C')
    fig.savefig(out/'rmsd.png',dpi=170);plt.close(fig)
    ranked=sorted(contacts,key=lambda r:-r['contact_occupancy'])[:20][::-1]
    fig,ax=plt.subplots(figsize=(9,7),layout='constrained')
    ax.barh([r['residue'] for r in ranked],[r['contact_occupancy'] for r in ranked],color='#387b92')
    ax.set_xlim(0,1);ax.set_xlabel(f'Fraction of observed frames with any heavy contact <= {args.contact_A:g} A')
    ax.set_title(f'{scope} | {args.case}\nTop residue contacts (not hydrogen-bond or affinity scores)')
    fig.savefig(out/'contact_occupancy.png',dpi=170);plt.close(fig)
    ranked=[r for r in hbonds if r['occupancy']>0][:15][::-1]
    fig,ax=plt.subplots(figsize=(11,max(3,.38*len(ranked)+1.6)),layout='constrained')
    if ranked:
        names=[f'{r["donor"].split("[")[0]} -> {r["acceptor"].split("[")[0]}' for r in ranked]
        ax.barh(names,[r['occupancy'] for r in ranked],color='#577d44');ax.set_xlim(0,1);ax.set_xlabel('Fraction of observed frames with qualified direct N/O H-bond')
    else:ax.text(.5,.5,'No qualified direct N/O hydrogen bonds in retained frames',ha='center',va='center');ax.set_axis_off()
    ax.set_title(f'{scope} | D-A <= {args.hbond_DA_A:g} A, H-A <= {args.hbond_HA_A:g} A, D-H-A >= {args.hbond_angle_deg:g} deg')
    fig.savefig(out/'hydrogen_bond_occupancy.png',dpi=170);plt.close(fig)
    limitations=details['warnings']+[
        'Occupancy is a fraction of retained sampled frames. Frames are correlated; no independent-trial confidence intervals or affinity conversion are claimed.',
        'Protein N/CA/C fitted to prepared input. Ligand RMSD uses fixed atom identities without symmetry correction; flexible loops can affect the global fit.',
        'Explicit donor H atoms and actual histidine protonation are used. Ligand N/O roles use RDKit BaseFeatures; protein N/O roles use residue chemistry and explicit bonded H.',
        'Direct N/O hydrogen bonds only. Water bridges, S-containing H-bonds, aromatic interactions, desolvation and binding free energies are outside this minimum analysis.',
        'Near H-bond CSV contains D/H/A triplets only in frames with D-A <=5 A; occupancy denominators still include all retained frames.',
        'One replicate or short pilot cannot establish convergence, improved affinity, pharmacological activity or dual-target function.'
    ]
    summary={'schema':1,'created_utc':datetime.now(timezone.utc).isoformat(),'case_id':args.case,'scope':scope,
             'software':{'python':sys.version,'numpy':np.__version__,'openmm':mm.__version__,'mdtraj':md.__version__,'matplotlib':matplotlib.__version__,'rdkit':rdBase.rdkitVersion},
             'coverage':{k:v for k,v in details.items() if k!='snapshot'},'input_segments':input_records,
             'rules':{'contact_distance_A':args.contact_A,'hbond_DA_A':args.hbond_DA_A,'hbond_HA_A':args.hbond_HA_A,'hbond_DHA_min_degrees':args.hbond_angle_deg,
                      'hbond_occupancy_unit':'donor-acceptor pair; any attached donor H meeting all three criteria in same frame','fit_atoms':len(backbone),'ligand_heavy_atoms':len(lheavy)},
             'ligand_chemical_roles':{'donors':[label(all_atoms[i]) for i in sorted(donors)],'acceptors':[label(all_atoms[i]) for i in sorted(acceptors)]},
             'rmsd':{key:{'mean':float(np.mean([r[key] for r in rmsd_rows])),'max':float(np.max([r[key] for r in rmsd_rows])),'final':rmsd_rows[-1][key]} for key in ('protein_backbone_RMSD_A','ligand_heavy_RMSD_A')},
             'qualified_hbond_pairs_with_nonzero_occupancy':sum(r['occupancy']>0 for r in hbonds),
             'limitations':limitations,'analysis_script_sha256':sha(Path(__file__))}
    (out/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (out/'source_manifest_snapshot.json').write_text(json.dumps(details['snapshot'],ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(out),'scope':scope,'retained_frames':done,'covered_sampling_ns':details['covered_sampling_ns'],
                      'rmsd':summary['rmsd'],'qualified_hbond_pairs':summary['qualified_hbond_pairs_with_nonzero_occupancy']},indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle-root',required=True,type=Path);p.add_argument('--case',required=True)
    p.add_argument('--run-dir',required=True,type=Path);p.add_argument('--out-dir',required=True,type=Path)
    p.add_argument('--legacy-engineering',action='store_true');p.add_argument('--allow-incomplete',action='store_true')
    p.add_argument('--chunk-frames',type=int,default=25)
    p.add_argument('--contact-A',type=float,default=4.0)
    p.add_argument('--hbond-DA-A',type=float,default=3.5);p.add_argument('--hbond-HA-A',type=float,default=2.5)
    p.add_argument('--hbond-angle-deg',type=float,default=135.0)
    p.add_argument('--highlight-residue',action='append',default=[],help='Additional protein residue ID to include in detailed tables; may be repeated. IDs are topology residue IDs, not canonical sequence indices.')
    args=p.parse_args()
    if args.chunk_frames<1 or min(args.contact_A,args.hbond_DA_A,args.hbond_HA_A)<=0 or not 0<args.hbond_angle_deg<=180:p.error('Invalid analysis parameters')
    analyze(args)


if __name__=='__main__':main()

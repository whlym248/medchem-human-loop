"""Read-only checks of an OpenMM engineering pilot (no Context or GPU needed).

Only Python stdlib, NumPy and OpenMM are required. This script never creates a
review approval. Geometry is checked at the final snapshot, not over the DCD.
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
import openmm as mm
from openmm import app, unit


PROTEIN = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP ASH GLH LYN CYX CYM HYP'.split())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def json_safe(value):
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key,item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def statistics(values):
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {'n_finite': 0}
    tail = a[len(a)//2:]
    return dict(n_finite=len(a), minimum=float(a.min()), maximum=float(a.max()),
                mean=float(a.mean()), sd=float(a.std()), final=float(a[-1]),
                second_half_mean=float(tail.mean()))


def csv_summary(path):
    if not path.is_file():
        return {'exists': False, 'path': str(path)}
    with path.open(encoding='utf-8-sig', newline='') as stream:
        first = stream.readline().lstrip('#')
        keys = next(csv.reader([first]))
        rows = list(csv.reader(stream))
    values = {key: [] for key in keys}
    bad = []
    placeholders = []
    for row_number, row in enumerate(rows, 2):
        if len(row) != len(keys):
            bad.append({'row': row_number, 'problem': 'column_count', 'count': len(row)})
            continue
        for key, raw in zip(keys, row):
            if raw.strip() in ('', '--') and key.startswith('Speed'):
                placeholders.append({'row': row_number, 'column': key, 'value': raw})
                continue
            try:
                number = float(raw)
            except ValueError:
                bad.append({'row': row_number, 'column': key, 'problem': 'not_numeric', 'value': raw})
                continue
            if not np.isfinite(number):
                bad.append({'row': row_number, 'column': key, 'problem': 'nonfinite', 'value': raw})
            values[key].append(number)
    steps = values.get('Step', [])
    return {'exists': True, 'path': str(path), 'sha256': sha(path), 'rows': len(rows),
            'columns': {key: statistics(v) for key, v in values.items()},
            'invalid_or_nonfinite': bad, 'expected_speed_placeholders': placeholders,
            'step_strictly_increasing': bool(len(steps) > 0 and np.all(np.diff(steps) > 0)),
            'scope': 'This attempt only; resumed trajectories must be trimmed by segment_manifest.json.'}


class MinimumImage:
    def __init__(self, box):
        self.box = np.asarray(box, dtype=float)
        self.inv = np.linalg.inv(self.box)
        self.orthogonal = bool(np.allclose(self.box @ self.box.T,
                                np.diag(np.diag(self.box @ self.box.T)), atol=1e-9))
        self.shifts = np.asarray(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)

    def __call__(self, delta):
        fraction = np.asarray(delta) @ self.inv
        centered = fraction - np.rint(fraction)
        base = centered @ self.box
        if self.orthogonal:
            return base
        # Nearest-image enumeration around the rounded image. OpenMM uses a
        # reduced periodic box; inspect 27 surrounding images for triclinic cells.
        best = base.copy()
        norm2 = np.sum(best*best, axis=-1)
        for shift in self.shifts:
            trial = (centered + shift) @ self.box
            test = np.sum(trial*trial, axis=-1)
            choose = test < norm2
            best = np.where(choose[..., None], trial, best)
            norm2 = np.minimum(norm2, test)
        return best


def atom_record(atom):
    return {'index_0': atom.index, 'chain': atom.residue.chain.id,
            'residue': atom.residue.name, 'resid': atom.residue.id,
            'atom': atom.name, 'element': atom.element.symbol if atom.element else None}


def bonds_and_constraints(system):
    bonds = {}
    for i in range(system.getNumConstraints()):
        a, b, length = system.getConstraintParameters(i)
        bonds[tuple(sorted((a,b)))] = float(length.value_in_unit(unit.nanometer))
    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            for i in range(force.getNumBonds()):
                a,b,length,k = force.getBondParameters(i)
                if float(k.value_in_unit(unit.kilojoule_per_mole/unit.nanometer**2)) > 0:
                    bonds[tuple(sorted((a,b)))] = float(length.value_in_unit(unit.nanometer))
    return bonds


def inspect_geometry(source, attempt, metadata):
    print('Reading input and final CIF (CPU only)...',file=sys.stderr,flush=True)
    initial = app.PDBxFile(str(source/'complex.cif'))
    final = app.PDBxFile(str(attempt/'final.cif'))
    atoms = list(initial.topology.atoms())
    final_atoms = list(final.topology.atoms())
    identity = [atom_record(a) for a in atoms] == [atom_record(a) for a in final_atoms]
    initial_xyz = np.load(source/'positions_nm.npy', allow_pickle=False)
    final_xyz = np.asarray(final.positions.value_in_unit(unit.nanometer), dtype=float)
    box = np.asarray(final.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer), dtype=float)
    report = {'atom_identity_order_matches_input': identity,
              'initial_coordinate_shape': list(initial_xyz.shape),
              'final_coordinate_shape': list(final_xyz.shape),
              'coordinates_all_finite': bool(np.isfinite(final_xyz).all()),
              'box_nm': box.tolist(), 'box_volume_nm3': float(np.linalg.det(box))}
    if not identity or final_xyz.shape != initial_xyz.shape or not np.isfinite(final_xyz).all() or np.linalg.det(box)<=0:
        return report
    mic = MinimumImage(box)
    report['minimum_image'] = 'orthogonal wrapping' if mic.orthogonal else '27 neighboring images of reduced triclinic box'
    ligand = set(metadata['ligand_atom_indices'])
    if any(i < 0 or i >= len(atoms) for i in ligand):
        raise ValueError('Ligand atom indices outside topology')
    heavy = np.array([a.index for a in atoms if a.element and a.element.atomic_number>1 and
                      (a.residue.name in PROTEIN or a.index in ligand)], dtype=int)
    heavy_set = set(heavy.tolist())
    protein_heavy = np.array([i for i in heavy if i not in ligand], dtype=int)
    ligand_heavy = np.array([i for i in heavy if i in ligand], dtype=int)
    report.update(solute_heavy_atoms=len(heavy), protein_heavy_atoms=len(protein_heavy), ligand_heavy_atoms=len(ligand_heavy))
    print('Reading force-field bond geometry...',file=sys.stderr,flush=True)
    system = mm.XmlSerializer.deserialize((source/'system.xml').read_text(encoding='utf-8'))
    if system.getNumParticles() != len(atoms):
        raise ValueError('System particles differ from topology atom count')
    bonds = bonds_and_constraints(system)
    adjacency = {i:set() for i in heavy}
    for a,b in bonds:
        if a in heavy_set and b in heavy_set:
            adjacency[a].add(b)
            adjacency[b].add(a)
    exclusions = {}
    for a in heavy:
        exclude = {a}|adjacency[a]
        for b in adjacency[a]:
            exclude.update(adjacency[b])
        exclusions[a] = exclude
    # All protein and ligand heavy atom pairs, excluding covalent 1-2 and 1-3.
    # Solvent and ions excluded deliberately: this is a solute geometry check.
    severe = []
    close_count = 0
    minimum = None
    pl_min = None
    print(f'Checking {len(heavy)} solute heavy atoms under PBC...',file=sys.stderr,flush=True)
    for start in range(0,len(heavy),128):
        ids = heavy[start:start+128]
        distance = np.linalg.norm(mic(final_xyz[ids,None,:]-final_xyz[heavy][None,:,:]),axis=-1)*10
        for row,a in enumerate(ids):
            allowed = np.array([int(b)>int(a) and b not in exclusions[a] for b in heavy])
            if not allowed.any():
                continue
            best = int(np.argmin(np.where(allowed,distance[row],np.inf)))
            record = (float(distance[row,best]),int(a),int(heavy[best]))
            if minimum is None or record[0]<minimum[0]:
                minimum=record
            close_count += int(np.sum(allowed & (distance[row]<1.8)))
            for col in np.flatnonzero(allowed & (distance[row]<1.5)):
                severe.append((float(distance[row,col]),int(a),int(heavy[col])))
        cross = np.array([[((int(a) in ligand)!=(int(b) in ligand)) for b in heavy] for a in ids])
        if cross.any():
            flat = int(np.argmin(np.where(cross,distance,np.inf)))
            row,col = np.unravel_index(flat,distance.shape)
            record=(float(distance[row,col]),int(ids[row]),int(heavy[col]))
            if pl_min is None or record[0]<pl_min[0]:
                pl_min=record
    def pair(r):
        return None if r is None else {'distance_A':r[0],'atom1':atom_record(atoms[r[1]]),'atom2':atom_record(atoms[r[2]])}
    report['nonbonded_solute_heavy_contacts']={
        'severe_threshold_A':1.5,'close_review_threshold_A':1.8,
        'severe_count':len(severe),'close_count_including_severe':close_count,
        'closest_nonbonded_solute_pair':pair(minimum),'closest_protein_ligand_pair':pair(pl_min),
        'severe_pairs_first_50':[pair(r) for r in sorted(severe)[:50]],
        'exclusions':'Exclude covalent 1-2 and 1-3; solvent/ions not examined. Thresholds flag gross overlap, not contact quality.'}
    deviations=[]
    broken=[]
    for (a,b),ideal in bonds.items():
        if a in heavy_set and b in heavy_set:
            actual=float(np.linalg.norm(mic(final_xyz[a]-final_xyz[b])))
            deviation=abs(actual-ideal)*10
            deviations.append(deviation)
            if deviation>0.30 or actual/ideal<0.65 or actual/ideal>1.5:
                broken.append({'atom1':atom_record(atoms[a]),'atom2':atom_record(atoms[b]),
                              'length_A':actual*10,'forcefield_equilibrium_A':ideal*10,'deviation_A':deviation})
    report['solute_heavy_bonds']={'n':len(deviations),'max_abs_deviation_from_forcefield_A':max(deviations,default=None),
                                 'review_count':len(broken),'review_bonds_first_50':broken[:50],
                                 'rule':'abs deviation >0.30 A or length/equilibrium outside [0.65,1.50]; flag for inspection, not covalent chemistry validation'}
    # Very short pilots permit a reference-based image mapping. This must not be
    # reused blindly for long trajectories whose atoms diffuse > half a box.
    mapped=initial_xyz+mic(final_xyz-initial_xyz)
    ca=np.array([i for i in protein_heavy if atoms[i].name=='CA'],dtype=int)
    fit=ca if len(ca)>=3 else protein_heavy
    initial_center=initial_xyz[fit].mean(axis=0)
    final_center=mapped[fit].mean(axis=0)
    x=mapped[fit]-final_center
    y=initial_xyz[fit]-initial_center
    print('Aligning final protein snapshot...',file=sys.stderr,flush=True)
    u,s,vt=np.linalg.svd(x.T@y)
    d=np.eye(3); d[2,2]=np.linalg.det(u@vt)
    rotation=u@d@vt
    aligned=(mapped-final_center)@rotation+initial_center
    rmsd=lambda ids:float(np.sqrt(np.mean(np.sum((aligned[ids]-initial_xyz[ids])**2,axis=1)))*10)
    report['final_snapshot_displacement']={'fit_atoms':'protein CA' if len(ca)>=3 else 'protein heavy',
        'protein_fit_rmsd_A':rmsd(fit),'protein_heavy_rmsd_A':rmsd(protein_heavy),
        'ligand_heavy_rmsd_after_protein_alignment_A':rmsd(ligand_heavy),
        'image_mapping':'Each final atom mapped to its nearest image relative to its input coordinate; suitable for a short engineering pilot only.',
        'interpretation':'Snapshot displacement is not pose stability, trajectory convergence or affinity.'}
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle-root',required=True,type=Path)
    p.add_argument('--case',required=True)
    p.add_argument('--attempt',required=True,type=Path,help='Directory containing final.cif and state.csv; v2 attempt_* or legacy seed directory')
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    root=args.bundle_root.resolve(); attempt=args.attempt.resolve()
    if '/' in args.case or '\\' in args.case or args.case in ('.','..'):
        p.error('Invalid case identifier')
    source=root/'systems'/args.case
    metadata=read_json(source/'metadata.json')
    if metadata.get('case_id')!=args.case:
        raise ValueError('Metadata case identity mismatch')
    report={'schema':1,'created_utc':datetime.now(timezone.utc).isoformat(),'case_id':args.case,
            'bundle_root':str(root),'attempt':str(attempt),'software':{'python':sys.version,'openmm':mm.__version__,'numpy':np.__version__},
            'automatic_approval_created':False,'scope':'Engineering quality screening only. Requires actual human/agent review. Final-snapshot geometry is not full trajectory analysis.',
            'critical_flags':[],'review_flags':[]}
    critical=report['critical_flags']; review=report['review_flags']
    report['input_sha256']={name:sha(source/name) for name in ('system.xml','complex.cif','positions_nm.npy','metadata.json')}
    if report['input_sha256']['system.xml']!=metadata['system_sha256']:
        critical.append('System XML hash does not match preparation metadata')
    status_file=next((attempt/name for name in ('attempt_status.json','run_status.json') if (attempt/name).is_file()),None)
    status=read_json(status_file) if status_file else {}
    report['run_status']=status
    if status.get('status')!='completed':
        critical.append('Run status missing or not completed')
    if status.get('case_id')!=args.case:
        critical.append('Run status case mismatch')
    final_energy=status.get('final_potential_energy_kj_mol')
    if not isinstance(final_energy,(int,float)) or not np.isfinite(final_energy):
        critical.append('Final potential energy missing/nonfinite in completed status')
    if 'system_sha256' in status and status['system_sha256']!=report['input_sha256']['system.xml']:
        critical.append('Legacy run system hash mismatch')
    identity_path=attempt.parent.parent/'run_identity.json'
    if identity_path.is_file():
        identity=read_json(identity_path)
        report['v2_run_identity_sha256']=sha(identity_path)
        if identity.get('input_sha256')!=report['input_sha256']:
            critical.append('v2 run identity input hashes differ from supplied bundle')
        canonical=json.dumps(identity,sort_keys=True,separators=(',',':'),ensure_ascii=False)
        identity_digest=hashlib.sha256(canonical.encode('utf-8')).hexdigest()
        if status.get('run_identity_sha256')!=identity_digest:
            critical.append('v2 attempt status is not bound to this run identity')
        if status.get('completed_absolute_step')!=identity.get('target_absolute_step'):
            critical.append('v2 completed absolute step differs from declared target')
    else:
        report['v2_run_identity_sha256']=None
        review.append('Legacy output: no v2 run identity is present; no new-run validation is implied')
    for phase,name in (('equilibration','equilibration.csv'),('sampling','state.csv')):
        summary=csv_summary(attempt/name)
        report[phase]=summary
        if not summary['exists']:
            review.append(f'{name} absent (resumed attempts may omit equilibration)')
            if phase=='sampling':critical.append('Sampling CSV absent')
            continue
        if summary['invalid_or_nonfinite']:
            critical.append(f'{name} has nonfinite/non-numeric fields')
        if not summary['rows'] or not summary['step_strictly_increasing']:
            critical.append(f'{name} empty or steps not strictly increasing')
        for required in ('Potential Energy (kJ/mole)','Temperature (K)','Density (g/mL)'):
            if summary['columns'].get(required,{}).get('n_finite',0)!=summary['rows']:
                critical.append(f'{name}: missing/invalid {required}')
    sampling=report['sampling'].get('columns',{})
    t=sampling.get('Temperature (K)',{})
    target=float(metadata.get('temperature_K',310))
    if 'mean' in t and abs(t['second_half_mean']-target)>15:
        review.append('Sampling second-half mean temperature differs from target by >15 K')
    density=sampling.get('Density (g/mL)',{})
    if 'mean' in density and not 0.7<=density['second_half_mean']<=1.3:
        review.append('Sampling second-half mean density outside broad 0.7–1.3 g/mL engineering screen')
    report['throughput']={k:status.get(k) for k in ('timed_steps','elapsed_s','ns_per_day')}
    report['throughput']['csv_speed']=sampling.get('Speed (ns/day)',{})
    report['throughput']['limitation']='Short timed sampling excludes environment setup, kernel compilation, minimization and initial equilibration; do not promise full campaign cost from this alone.'
    if not isinstance(status.get('ns_per_day'),(int,float)) or not np.isfinite(status.get('ns_per_day',float('nan'))) or status.get('ns_per_day',0)<=0:
        review.append('No finite positive status ns/day')
    final_path=attempt/'final.cif'
    if final_path.is_file():
        report['final_cif_sha256']=sha(final_path)
        geometry=inspect_geometry(source,attempt,metadata)
        report['geometry']=geometry
        if not geometry['atom_identity_order_matches_input'] or not geometry['coordinates_all_finite'] or geometry['box_volume_nm3']<=0 or geometry['initial_coordinate_shape']!=geometry['final_coordinate_shape']:
            critical.append('Invalid final geometry identity/coordinates/box')
        if geometry.get('nonbonded_solute_heavy_contacts',{}).get('severe_count',0)>0:
            critical.append('Final snapshot has solute heavy nonbonded overlaps <1.5 A')
        if geometry.get('nonbonded_solute_heavy_contacts',{}).get('close_count_including_severe',0)>0:
            review.append('Inspect final solute heavy nonbonded contacts <1.8 A')
        if geometry.get('solute_heavy_bonds',{}).get('review_count',0)>0:
            review.append('Inspect flagged final solute heavy bond lengths')
    else:
        critical.append('final.cif missing')
    report['result']='critical_flags_require_resolution' if critical else 'screening_completed_requires_review'
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(json_safe(report),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(json_safe({'result':report['result'],'out':str(args.out.resolve()),'critical_flags':critical,
                      'review_flags':review,'ns_per_day':status.get('ns_per_day'),
                      'final_geometry':report.get('geometry',{})}),ensure_ascii=False,indent=2,allow_nan=False))
    return 2 if critical else 0


if __name__=='__main__':
    try:
        sys.exit(main())
    except Exception as error:
        print(f'Quality-check script failed: {type(error).__name__}: {error}',file=sys.stderr)
        sys.exit(3)

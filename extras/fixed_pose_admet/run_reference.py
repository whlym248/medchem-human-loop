"""Reference interface: score two supplied poses and run pretrained ADMET-AI.

Requires an explicitly prepared external input directory and existing tools and
weights. No MD, training, pose search or downloads are performed by this script.
See README.md for the input contract and validation limits.
"""
import os
os.environ['CUDA_VISIBLE_DEVICES']='-1'
for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:
    os.environ[name]='2'
from pathlib import Path
from datetime import datetime, timezone
import argparse, csv, hashlib, importlib.metadata, json, math, platform, re, shutil, subprocess, sys, time

ROOT=Path(__file__).resolve().parent
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def require(condition, message):
    if not condition:raise ValueError(message)

def relative_path(value):
    p=(ROOT/value).resolve()
    if not p.is_relative_to(ROOT):raise ValueError(f'Input outside example directory: {value}')
    return p

def main():
    global ROOT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True,help='Prepared input directory; never supplied in the public repository')
    parser.add_argument('--config',default='config.local.json')
    parser.add_argument('--output',default='runs/'+datetime.now().strftime('%Y%m%d_%H%M%S'))
    args=parser.parse_args()
    ROOT=args.root.resolve()
    if not ROOT.is_dir():raise ValueError('Prepared input directory must exist')
    config_path=Path(args.config);config_path=config_path if config_path.is_absolute() else ROOT/config_path
    config=load(config_path);protocol=load(ROOT/'protocol.json');manifest=load(ROOT/'input_manifest.json')
    records=load(ROOT/'data/candidates.json')
    require(len(records)==2, 'This reference contract requires exactly two candidates')
    ids=[r['candidate_id'] for r in records]
    require(len(set(ids))==2 and all(re.fullmatch(r'[A-Za-z0-9_-]+', x) for x in ids), 'Candidate IDs must be distinct safe filename stems')
    require(protocol['scoring_modes']==['vina','vinardo'], 'Expected Vina and Vinardo score-only modes')
    require(protocol['cpu_per_scoring_job']==1, 'Reference scoring uses one CPU worker')
    require(bool(protocol['expected_versions']), 'Record and enforce dependency versions')
    for key in ['center','size']:
        require(len(protocol['box'][key])==3 and all(isinstance(v,(int,float)) and math.isfinite(v) for v in protocol['box'][key]), 'Box must have three finite values')
    require(all(v>0 for v in protocol['box']['size']), 'Box sizes must be positive')
    required_files={'data/candidates.json','protocol.json','expected_models.json'}
    for record in records:
        required_files.update(record[k] for k in ['pose_pdbqt','structure_file','receptor_pdbqt','receptor_structure_file'])
    require(required_files <= {item['path'] for item in manifest['files']}, 'Input manifest must cover configuration and all referenced structures')
    out=Path(args.output);out=out if out.is_absolute() else ROOT/out
    out.mkdir(parents=True,exist_ok=False)
    logs=out/'logs';logs.mkdir()
    started=datetime.now(timezone.utc).isoformat();t0=time.time()
    save(out/'status.json',{'status':'running','started_utc':started})
    try:
        for item in manifest['files']:
            p=relative_path(item['path'])
            if sha(p)!=item['sha256']:raise ValueError(f'Input hash mismatch: {item["path"]}')
        requested=config['vina_executable'];vina=shutil.which(requested) or requested
        if not Path(vina).is_file():raise FileNotFoundError('Configure an existing Vina executable; this script does not download it.')
        if config.get('vina_expected_sha256') and sha(vina)!=config['vina_expected_sha256']:raise ValueError('Vina hash mismatch')
        v=subprocess.run([vina,'--version'],capture_output=True,text=True,timeout=30,check=True)
        (logs/'vina_version.log').write_text(v.stdout+v.stderr,encoding='utf-8')
        metadata={'started_utc':started,'python':sys.version,'python_executable':sys.executable,'os':platform.platform(),
            'operation':'Fresh score_only evaluations and fresh pretrained ADMET inference; no previous output reuse',
            'vina_executable':str(vina),'vina_sha256':sha(vina),'vina_version':v.stdout.strip(),
            'input_files':manifest['files'],'protocol_sha256':sha(ROOT/'protocol.json'),'main_py_sha256':sha(Path(__file__)),
            'config_sha256':sha(config_path),'seed':protocol['seed'],'device':'cpu','models_downloaded':False,'trained_or_finetuned':False}
        score_records=[]
        for r in records:
            for method in protocol['scoring_modes']:
                command=[vina,'--receptor',str(relative_path(r['receptor_pdbqt'])),'--ligand',str(relative_path(r['pose_pdbqt'])),
                    '--scoring',method,'--score_only','--cpu',str(protocol['cpu_per_scoring_job']),'--seed',str(protocol['seed'])]
                for key in ['center','size']:
                    for axis,value in zip('xyz',protocol['box'][key]):command.extend([f'--{key}_{axis}',str(value)])
                call_start=time.time();run=subprocess.run(command,capture_output=True,text=True,timeout=120)
                text=run.stdout+'\n'+run.stderr
                logfile=logs/f'{r["candidate_id"]}_{method}_score_only.log';logfile.write_text(text,encoding='utf-8')
                match=re.search(r'Estimated Free Energy of Binding\s*:\s*([-\d.]+)',text) or re.search(r'Affinity:\s*([-\d.]+)',text)
                if run.returncode or not match:raise RuntimeError(f'Scoring failed: {logfile}')
                score=float(match[1]);require(math.isfinite(score), 'Nonfinite score')
                score_records.append({'candidate_id':r['candidate_id'],'method':method,'operation':'score_only','score_kcal_mol':score,
                    'command':command,'exit_code':run.returncode,'elapsed_seconds':time.time()-call_start,'log_file':str(logfile.relative_to(out)),
                    'log_sha256':sha(logfile)})
                print(f'{r["candidate_id"]} {method} fresh score: {score}',flush=True)
        save(out/'fresh_scoring.json',score_records)

        # Existing local weights are hash-checked before the library loads them.
        import numpy as np
        import torch
        torch.set_num_threads(2);torch.set_num_interop_threads(1);torch.manual_seed(protocol['seed']);np.random.seed(protocol['seed'])
        from rdkit import Chem
        from admet_ai import ADMETModel
        from admet_ai.constants import DEFAULT_MODELS_DIR
        modeldir=Path(config.get('admet_models_dir') or DEFAULT_MODELS_DIR)
        expected=load(ROOT/'expected_models.json');weight_records=[]
        for item in expected['files']:
            path=(modeldir/item['relative_path']).resolve()
            require(path.is_relative_to(modeldir.resolve()), 'Model path escapes model directory')
            if not path.is_file():raise FileNotFoundError(f'Local model missing; no download attempted: {path}')
            if sha(path)!=item['sha256']:raise ValueError(f'Weight hash mismatch: {path}')
            weight_records.append({**item,'loaded_path':str(path)})
        actual_weights=sorted(str(p.relative_to(modeldir)).replace('\\','/') for p in modeldir.rglob('*.pt'))
        require(bool(expected['files']) and actual_weights==sorted(x['relative_path'] for x in expected['files']), 'Model manifest must cover all installed .pt weights')
        versions={name:importlib.metadata.version(name) for name in ['admet-ai','chemprop','torch','rdkit','numpy','pandas','lightning']}
        for name,value in protocol['expected_versions'].items():
            if versions[name]!=value:raise ValueError(f'Expected {name} {value}, got {versions[name]}')
        smiles=[]
        for r in records:
            mol=Chem.MolFromSmiles(r['smiles']);require(mol is not None, 'Invalid SMILES')
            require(Chem.MolToSmiles(mol,isomericSmiles=True)==r['canonical_isomeric_smiles'], 'Canonical isomeric SMILES mismatch')
            require(Chem.GetFormalCharge(mol)==r['formal_charge'], 'Formal charge mismatch');smiles.append(r['smiles'])
        model=ADMETModel(models_dir=modeldir,num_workers=0,drugbank_path=None)
        require(model.device=='cpu', 'Expected CPU inference')
        require([len(x) for x in model.model_lists]==[5,5], 'Expected two five-model ensembles')
        print('Running fresh ADMET-AI inference for two inputs on CPU',flush=True)
        pred=model.predict(smiles)
        require(list(pred.index)==smiles, 'Prediction row order mismatch')
        pred.insert(0,'candidate_id',[r['candidate_id'] for r in records]);pred.insert(1,'input_smiles',smiles)
        pred.to_csv(out/'fresh_admet_full.csv',index=False,encoding='utf-8')
        score_lookup={(r['candidate_id'],r['method']):r['score_kcal_mol'] for r in score_records}
        results=[]
        for i,r in enumerate(records):
            p=pred.iloc[i]
            values={key:float(p[key]) for key in ['molecular_weight','tpsa','logP','Solubility_AqSolDB','hERG','Caco2_Wang']}
            require(all(math.isfinite(x) for x in values.values()), 'Nonfinite prediction')
            results.append({'candidate_id':r['candidate_id'],'record_group':r.get('record_group','user_supplied'),'target':r['target'],'smiles':r['smiles'],
                'formal_charge':r['formal_charge'],'structure_file':r['structure_file'],'receptor_structure_file':r['receptor_structure_file'],
                'source_docking_seed':r['source_docking_seed'],'source_pose_rank':r['source_pose_rank'],
                'vina_score_only_kcal_mol':score_lookup[(r['candidate_id'],'vina')],
                'vinardo_score_only_kcal_mol':score_lookup[(r['candidate_id'],'vinardo')],
                'molecular_weight_Da':values['molecular_weight'],'tpsa_A2':values['tpsa'],'clogp_descriptor':values['logP'],
                'predicted_logS_AqSolDB':values['Solubility_AqSolDB'],'hERG_model_score':values['hERG'],'Caco2_Wang_model_output':values['Caco2_Wang'],
                'scoring_software_version':v.stdout.strip(),'ai_model':'ADMET-AI','ai_version':versions['admet-ai'],'chemprop_version':versions['chemprop'],
                'notes':'Fresh fixed-pose rescore and pretrained inference; not affinity, activity, MD, new design, or patient-risk validation.'})
        with (out/'results.csv').open('w',encoding='utf-8',newline='') as fh:
            writer=csv.DictWriter(fh,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
        metadata.update(actual_package_versions=versions,model_weights=weight_records,model_ensemble_sizes=[len(x) for x in model.model_lists],
            admet_drugbank_percentiles=False,completed_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-t0,
            fresh_scoring_evaluations=len(score_records),fresh_admet_inputs=len(pred),final_results_sha256=sha(out/'results.csv'))
        save(out/'runtime_metadata.json',metadata)
        save(out/'status.json',{'status':'complete','started_utc':started,'completed_utc':metadata['completed_utc'],
            'elapsed_seconds':metadata['elapsed_seconds'],'fresh_score_evaluations':len(score_records),'fresh_admet_inputs':len(pred)})
        print(f'COMPLETE {out}',flush=True)
    except Exception as e:
        save(out/'status.json',{'status':'failed','started_utc':started,'error_type':type(e).__name__,'error':str(e)})
        raise

if __name__=='__main__':main()

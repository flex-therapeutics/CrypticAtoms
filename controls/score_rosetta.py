from pathlib import Path
import os,sys,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from concurrent.futures import ProcessPoolExecutor
ROOT=Path(__file__).resolve().parent;REVIEW=ROOT.parent;REPO=REVIEW;sys.path.insert(0,str(REPO))
from benchmark.metrics import load_protein_chain,superimpose_on_pocket_ca,all_atom_rmsd_over_residues
T={t['name']:t for t in json.loads((REPO/'benchmark/targets.json').read_text())};manifest=json.loads((ROOT/'inputs.json').read_text())
REF={}
for name,t in T.items():
 a=load_protein_chain(REPO/t['apo']['path'],t['apo']['chain']);h=load_protein_chain(REPO/t['holo']['path'],t['holo']['chain']);REF[name]=(a,h)
def score(job):
 arm,name,key,i=job;t=T[name];out={'plot_method':'rosetta_'+arm,'target':name,'sample_idx':i,'status':'ok','error':''}
 try:
  path=ROOT/'outputs'/arm/key/f'sample_{i:03d}.pdb';p=load_protein_chain(path,t['apo']['chain']);a,h=REF[name];pocket=t['pocket']['residues_5A'];moving=t['pocket']['moving_residues_2A'];pa,ca_a,na=superimpose_on_pocket_ca(p,a,pocket);ph,ca_h,nh=superimpose_on_pocket_ca(p,h,pocket)
  da=all_atom_rmsd_over_residues(pa,a,moving);dh=all_atom_rmsd_over_residues(ph,h,moving);D=t['baselines']['align_pocket']['aa_rmsd_moving']
  assert np.isfinite(da) and np.isfinite(dh)
  out.update(aa_rmsd_moving_vs_apo=da,aa_rmsd_moving_vs_holo=dh,D=D,tier1=bool(dh<D and dh<da),pocket_ca_rmsd_vs_apo=ca_a,pocket_ca_rmsd_vs_holo=ca_h,n_align_apo=na,n_align_holo=nh,structure_path=str(path))
 except Exception as e:out.update(status='error',error=str(e))
 return out
if __name__=='__main__':
 jobs=[(arm,t['target'],t['input_key'],i) for t in manifest for arm in ['repack','fastrelax'] for i in range(30) if (ROOT/'outputs'/arm/t['input_key']/f'sample_{i:03d}.json').exists()]
 with ProcessPoolExecutor(max_workers=16) as pool:out=list(pool.map(score,jobs))
 if not out:raise SystemExit('No complete samples')
 raw=pd.DataFrame(out);raw.to_parquet(ROOT/'selected_structural_scores.parquet',index=False);rows=[]
 for t in manifest:
  for arm in ['repack','fastrelax']:
   name=t['target'];slug='rosetta_'+arm;q=raw[(raw.target==name)&(raw.plot_method==slug)];v=q[q.status=='ok'];n=len(v);D=T[name]['baselines']['align_pocket']['aa_rmsd_moving']
   if n:
    da=v.aa_rmsd_moving_vs_apo.to_numpy();dh=v.aa_rmsd_moving_vs_holo.to_numpy();base=(dh<D)&(dh<da);x=(D*D+da*da-dh*dh)/(2*D);y2=da*da-x*x;assert y2.min()>-1e-4,(name,y2.min());y=np.sqrt(np.maximum(y2,0))
   rows.append(dict(method=slug,target=name,D=D,n_saved=len(q),n_selected=len(q),n_valid=n,tier1=int(base.any()) if n else 0,tier2=int((base&(dh<1.5)).any()) if n else 0,tier3=int((base&(dh<1.)).any()) if n else 0,n_tier1=int(base.sum()) if n else 0,map_rg=float(np.sqrt(np.mean((x-x.mean())**2+(y-y.mean())**2))) if n else np.nan,sigma_x=float(np.std(x,ddof=1)) if n>1 else np.nan,sigma_y=float(np.std(y,ddof=1)) if n>1 else np.nan))
 summary=pd.DataFrame(rows);summary.to_csv(ROOT/'structural_per_entry.csv',index=False);print(raw.groupby(['plot_method','status']).size());print(summary.groupby('method')[['tier1','tier2','tier3','n_valid']].sum());print(raw[raw.status!='ok'][['target','error']].drop_duplicates().to_string(index=False))

from pathlib import Path
import sys,json,numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment
import biotite.structure.io as io
r=Path(__file__).parent;w=r.parent;sys.path.insert(0,str(w/'validity'));import compute as V
from benchmark.metrics import superimpose_on_pocket_ca,all_atom_rmsd_over_residues,_kabsch_align_residue,match_atoms
methods=['rosetta_repack','rosetta_fastrelax','prody_rfd3','bioemu_refine','boltz2','esmfold2'];targets=['cAbl','HIV_RT__3V81','GLTP'];s=pd.concat([pd.read_parquet(w/'analysis/selected_structural_scores.parquet'),pd.read_parquet(w/'controls/selected_structural_scores.parquet')],ignore_index=True);s=s[s.target.isin(targets)&s.plot_method.isin(methods)].copy();manifest=pd.read_csv(w/'validity/input_manifest.csv');rows=[];details=[];refs=[];cache={};refarrays={};bb={'N','CA','C','O'}
def split(a,b,moving):
 vals=[]
 for rid in moving:
  aa=a[(a.res_id==rid)&(a.atom_name!='OXT')];hh=b[(b.res_id==rid)&(b.atom_name!='OXT')];assert set(aa.atom_name)==set(hh.atom_name),(rid,'incomplete atoms');assert set(aa.res_name)==set(hh.res_name),(rid,'residue mismatch')
  ab,hb=match_atoms(aa[np.isin(aa.atom_name,list(bb))],hh[np.isin(hh.atom_name,list(bb))]);bsq=((ab.coord-hb.coord)**2).sum(1);al=_kabsch_align_residue(aa,hh);sc=[]
  for el in set(aa.element):
   ai=np.where((aa.element==el)&~np.isin(aa.atom_name,list(bb)))[0];hi=np.where((hh.element==el)&~np.isin(hh.atom_name,list(bb)))[0];assert len(ai)==len(hi)
   if len(ai):
    ii,jj=linear_sum_assignment(np.linalg.norm(al.coord[ai,None,:]-hh.coord[hi][None,:,:],axis=-1));sc.extend(((aa.coord[ai[ii]]-hh.coord[hi[jj]])**2).sum(1).tolist())
  vals.append(dict(res_id=int(rid),res_name=str(hh.res_name[0]),bb_sum=float(bsq.sum()),bb_n=len(bsq),sc_sum=float(sum(sc)),sc_n=len(sc)))
 return vals
for target in targets:
 t=V.T[target];a=V.load_protein_chain(V.REPO/t['apo']['path'],t['apo']['chain']);h=V.load_protein_chain(V.REPO/t['holo']['path'],t['holo']['chain']);p=t['pocket']['residues_5A'];mov=t['pocket']['moving_residues_2A'];ah,_,_=superimpose_on_pocket_ca(a,h,p);refarrays[target]=(ah,h);splitmov=[rid for rid in mov if set(a[a.res_id==rid].atom_name)==set(h[h.res_id==rid].atom_name) and set(a[a.res_id==rid].res_name)==set(h[h.res_id==rid].res_name)];rr=split(ah,h,splitmov);refs.append(dict(target=target,n_moving=len(mov),n_diagnostic=len(splitmov),excluded_residues=','.join(map(str,set(mov)-set(splitmov))),D=all_atom_rmsd_over_residues(ah,h,mov),bb_rmsd=np.sqrt(sum(z['bb_sum'] for z in rr)/sum(z['bb_n'] for z in rr)),sc_rmsd=np.sqrt(sum(z['sc_sum'] for z in rr)/sum(z['sc_n'] for z in rr))))
 bench='cp19' if 'reference_bundle/cp19/' in t['holo']['path'] else 'cp15';bundle=V.load(V.REPO/f'docking/reference_bundle/{bench}/{target}/holo.pdb');ids=list(dict.fromkeys(bundle[bundle.atom_name=='CA'].res_id.tolist()))
 for q in s[s.target.eq(target)].to_dict('records'):
  m=q['plot_method'];idx=int(q['sample_idx']);paths=[str(q.get(k,'')) for k in ['structure_path','struct_path']];found=manifest[(manifest.method==m)&(manifest.target==target)&(manifest.sample_idx==idx)];paths += found.local_path.tolist();path=next((Path(v) for v in paths if Path(v).is_file()),None)
  if m=='bbflow_refine':
   path=w/'bbflow_refresh_20260916/scored_pdbs'/target/f'sample_{idx}.pdb'
  assert path,(m,target,idx)
  x=V.load(path)
  if (m,target) in V.VERIFIED_CASES:x,_=V.remap_structure(x,m,target,sample_idx=idx,source_path=str(path),reference_chain=t['apo']['chain'],allow_substitutions=True)
  elif m in ['bioemu_refine','boltz2','esmfold2']:x.res_id=np.array([ids[int(v)-1] if 1<=int(v)<=len(ids) else -999999 for v in x.res_id]);x=x[x.res_id!=-999999]
  xh,_,n=superimpose_on_pocket_ca(x,h,p);xa,_,_=superimpose_on_pocket_ca(x,a,p);dh=all_atom_rmsd_over_residues(xh,h,mov);da=all_atom_rmsd_over_residues(xa,a,mov);stored=float(q['aa_rmsd_moving_vs_holo']);assert abs(dh-stored)<.02,(m,target,idx,dh,stored);vals=split(xh,h,splitmov);cache[(target,m,idx)]=xh
  record=dict(target=target,method=m,sample_idx=idx,Dsh=dh,Dsa=da,D=refs[-1]['D'],tier1=bool(dh<refs[-1]['D'] and dh<da),bb_rmsd=np.sqrt(sum(z['bb_sum'] for z in vals)/sum(z['bb_n'] for z in vals)),sc_rmsd=np.sqrt(sum(z['sc_sum'] for z in vals)/sum(z['sc_n'] for z in vals)),structure_path=str(path));rows.append(record)
  for z in vals:details.append(dict(target=target,method=m,sample_idx=idx,**z))
 print('DONE',target,flush=True)
d=pd.DataFrame(rows);d.to_csv(r/'per_sample.csv',index=False);pd.DataFrame(details).to_csv(r/'per_residue.csv',index=False);pd.DataFrame(refs).to_csv(r/'references.csv',index=False);selected=[];(r/'structures').mkdir(exist_ok=True)
for target in targets:
 for name,ar in zip(['apo','holo'],refarrays[target]):io.save_structure(str(r/'structures'/f'{target}_{name}.pdb'),ar)
 for group,ml in [('physics',methods[:3]),('learned',methods[3:])]:
  pool=d[d.target.eq(target)&d.method.isin(ml)];pool=pool[pool.tier1 & pool.method.eq({'cAbl':'boltz2','HIV_RT__3V81':'esmfold2','GLTP':'bioemu_refine'}[target])] if group=='learned' else pool;q=pool.sort_values(['Dsh','method','sample_idx']).iloc[0];rec=q.to_dict();rec['group']=group;selected.append(rec);io.save_structure(str(r/'structures'/f'{target}_{group}.pdb'),cache[(target,q.method,int(q.sample_idx))])
pd.DataFrame(selected).to_csv(r/'selected_examples.csv',index=False);print(pd.DataFrame(refs).to_string(index=False));print(pd.DataFrame(selected).to_string(index=False));print(d.groupby(['target','method']).tier1.sum().to_string())

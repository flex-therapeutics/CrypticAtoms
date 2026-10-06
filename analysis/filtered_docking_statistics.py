from pathlib import Path
import sys,itertools,pandas as pd
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R))
from review_docking_statistics import inference,holm,SEED
p=pd.read_csv(R/'clash_filtered_docking_per_target.csv');assert p.n_pending.sum()==0
w=p.pivot(index='target',columns='method',values='delta_apo')
generators=['boltz2','esmfold2','bioemu_refine','bbflow_refine','confornet_nomsa','rfd3_fix7','apo2mol','molgen_ns1.0'];controls=['rosetta_repack','rosetta_fastrelax']
families=[('vs_apo',[(m,'apo') for m in generators]),('pairwise',list(itertools.combinations(generators,2))),('prody_supporting',[('prody_rfd3',m) for m in ['apo']+generators]),('rosetta_supporting',[(m,c) for m in generators+['prody_rfd3'] for c in controls]+[(c,'apo') for c in controls]+[tuple(controls)])]
out=[]
for j,(name,pairs) in enumerate(families):
 rows=[]
 for i,(a,b) in enumerate(pairs):
  v=w[a] if b=='apo' else w[a]-w[b]
  rows.append(dict(method=a,comparator=b,**inference(v,SEED+j*1000+i)))
 d=pd.DataFrame(rows);d['cluster_signflip_holm_p']=holm(d.cluster_signflip_p);d['wilcoxon_holm_p']=holm(d.wilcoxon_p);d['family']=name;d['n_family_tests']=len(d)
 d.to_csv(R/f'review_docking_statistics_{name}.csv',index=False);out.append(d)
 print(name,len(d),'significant:',d[d.cluster_signflip_holm_p.lt(.05)][['method','comparator','mean_difference','cluster_signflip_holm_p']].to_string(index=False))
pd.concat(out).to_csv(R/'review_docking_statistics_all_comparisons.csv',index=False)

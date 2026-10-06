"""Paired tier recovery with protein clustering and complete-case missingness."""
from pathlib import Path
from itertools import combinations
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import binomtest

A=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--include-controls',action='store_true')
parser.add_argument('--no-figure',action='store_true')
args=parser.parse_args()
s=pd.read_csv(A/'structural_per_entry.csv');meta=pd.read_csv(A/'reference_metadata.csv')
# Preserve the original planned family; a later classical arm is supporting.
PRIMARY_PIPELINES=['boltz2','esmfold2','bioemu_refine','bbflow_refine','confornet_nomsa','rfd3_fix7','apo2mol','molgen_ns1.0']
CLASSICAL_PIPELINES=['prody_rfd3']
methods=json.loads((A/'methods.json').read_text())
metadata={m[0]:m for m in methods}
assert len(metadata)==len(methods), 'Duplicate method metadata'
assert set(PRIMARY_PIPELINES)<=set(metadata), 'Original family members must remain present'
assert set(metadata)<=set(PRIMARY_PIPELINES+CLASSICAL_PIPELINES), 'Define inference families explicitly for new methods'
methods=[metadata[m] for m in PRIMARY_PIPELINES+CLASSICAL_PIPELINES if m in metadata]
generators=PRIMARY_PIPELINES  # Original eight; excludes the classical ProDy sampler.
classical=[m for m in CLASSICAL_PIPELINES if m in metadata]
sampling_pipelines=PRIMARY_PIPELINES+classical
for slug in classical:
 arm=s[s.method==slug]
 assert len(arm)==len(meta) and arm.target.is_unique and set(arm.target)==set(meta.target), 'Require the full classical-pipeline target table'
 assert arm.n_selected.eq(30).all(), 'Classical-pipeline sample selection is not complete'

control_slugs=['rosetta_repack','rosetta_fastrelax']
if args.include_controls:
 path=A.parent/'controls/structural_per_entry.csv'
 if not path.exists():raise SystemExit('Controls missing; outputs preserved.')
 controls=pd.read_csv(path)
 aliases={'repack':'rosetta_repack','rosetta_repack':'rosetta_repack','fastrelax':'rosetta_fastrelax','rosetta_fastrelax':'rosetta_fastrelax','rosettarelax':'rosetta_fastrelax'}
 controls['method']=controls.method.map(aliases)
 additions=[]
 for slug,label,color in zip(control_slugs,['Rosetta repack','Rosetta FastRelax'],['#929292','#4F4F4F']):
  arm=controls[controls.method==slug].copy()
  complete=(len(arm)==len(meta) and arm.target.is_unique and set(arm.target)==set(meta.target) and arm.n_valid.eq(30).all() and arm[['tier1','tier2','tier3']].isin([0,1]).all().all())
  if 'n_selected' in arm:complete=complete and arm.n_selected.eq(30).all()
  if not complete:raise SystemExit(f'{label}: require all 21 targets with 30 valid samples; outputs preserved.')
  if slug in set(s.method):raise SystemExit('Control already present in primary data.')
  additions.append(arm);methods.append([slug,label,slug,color])
 s=pd.concat([s,*additions],ignore_index=True)
slugs=[m[0] for m in methods];labels={m[0]:m[1] for m in methods};order=meta.target.tolist()
assert len(order)==21 and meta.target.str.split('__').str[0].nunique()==17
assert not s.duplicated(['method','target']).any()
# Absent/invalid predictions are missing, never unsuccessful recovered-pocket outcomes.
s.loc[s.n_valid.fillna(0).le(0),['tier1','tier2','tier3']]=np.nan
seed=732;bootstrap_cache={};sign_cache={}
def cluster_arrays(targets):
 proteins=[t.split('__')[0] for t in targets];groups=list(dict.fromkeys(proteins))
 G=np.array([[p==g for p in proteins] for g in groups],dtype=int)
 key=tuple(targets)
 if key not in bootstrap_cache:
  rng=np.random.default_rng(seed);bootstrap_cache[key]=rng.multinomial(len(groups),np.full(len(groups),1/len(groups)),size=20000)
 return G,bootstrap_cache[key]
def interval(values,targets):
 G,boot=cluster_arrays(targets);reps=(boot@(G@values))/(boot@G.sum(axis=1))
 return np.quantile(reps,[.025,.975]),G

def holm_fixed(pvalues):
 # A non-test occupies its planned family slot as p=1 internally, but remains NaN in output.
 p=np.array(pvalues,dtype=float);finite=np.isfinite(p);working=np.where(finite,p,1.)
 idx=np.argsort(working);adj=np.maximum.accumulate((len(p)-np.arange(len(p)))*working[idx]);out=np.empty(len(p));out[idx]=np.minimum(adj,1);out[~finite]=np.nan;return out
rows=[];summ=[]
for tier in range(1,4):
 y=s.pivot(index='target',columns='method',values=f'tier{tier}').reindex(index=order,columns=slugs)
 for slug in slugs:
  valid=y[slug].dropna();ci,G=interval(valid.to_numpy(),valid.index) if len(valid) else ([np.nan,np.nan],np.empty((0,0)))
  summ.append(dict(method=slug,tier=tier,n_recovered=int(valid.sum()),n_evaluable=len(valid),n_missing=21-len(valid),n_proteins=len(G),recovery_fraction=valid.mean(),cluster_bootstrap_ci_low=ci[0],cluster_bootstrap_ci_high=ci[1]))
 for a,b in combinations(slugs,2):
  pair=y[[a,b]].dropna();v=(pair[a]-pair[b]).to_numpy();n=len(v)
  wins=int((v>0).sum());losses=int((v<0).sum());discordant=wins+losses
  status='tested' if discordant else ('identical' if n else 'no_evaluable_pairs')
  ci,G=interval(v,pair.index) if n else ([np.nan,np.nan],np.empty((0,0)));g=len(G)
  cp=mp=np.nan
  if discordant:
   if g not in sign_cache:sign_cache[g]=2*((np.arange(2**g,dtype=np.uint32)[:,None]>>np.arange(g))&1).astype(np.int8)-1
   cp=float(np.mean(np.abs(sign_cache[g]@(G@v))>=abs(v.sum())))
   mp=binomtest(wins,discordant,.5).pvalue
  rows.append(dict(tier=tier,method_a=a,method_b=b,status=status,n_evaluable=n,n_missing=21-n,n_proteins=g,a_only=wins,b_only=losses,n_discordant=discordant,paired_difference=v.mean() if n else np.nan,cluster_bootstrap_ci_low=ci[0],cluster_bootstrap_ci_high=ci[1],cluster_signflip_p=cp,mcnemar_exact_p=mp))
r=pd.DataFrame(rows);n_tests=len(r)
for col in ['cluster_signflip_p','mcnemar_exact_p']:
 r[col+'_holm_all_pairs']=holm_fixed(r[col])
 r[col+'_holm_controls_within_tier']=np.nan
 r[col+'_holm_prody_supporting']=np.nan
if args.include_controls:
 for tier in range(1,4):
  mask=(r.tier==tier)&r.method_a.isin(generators)&r.method_b.isin(control_slugs)
  assert mask.sum()==16
  for col in ['cluster_signflip_p','mcnemar_exact_p']:
   r.loc[mask,col+'_holm_controls_within_tier']=holm_fixed(r.loc[mask,col])
 # Six later-added ProDy/control hypotheses across all three tiers form a
 # distinct supporting family; they never change the original 16-per-tier family.
 support=r.method_a.isin(classical)&r.method_b.isin(control_slugs)
 if classical:
  assert support.sum()==6
  for col in ['cluster_signflip_p','mcnemar_exact_p']:
   r.loc[support,col+'_holm_prody_supporting']=holm_fixed(r.loc[support,col])
r.to_csv(A/'review_statistics_pairwise.csv',index=False)
pd.DataFrame(summ).to_csv(A/'review_statistics_recovery.csv',index=False)
cmask=r.method_a.isin(generators)&r.method_b.isin(control_slugs)
r.loc[cmask].to_csv(A/'review_statistics_control_comparisons.csv',index=False)
prody_mask=r.method_a.isin(classical)&r.method_b.isin(control_slugs)
r.loc[prody_mask].to_csv(A/'review_statistics_prody_control_comparisons.csv',index=False)
# One compact Tier-1 table: eight generators against each of two controls.
tex=[]
if args.include_controls:
 def pformat(p):
  return '--' if pd.isna(p) else (r'$<0.001$' if p<.001 else f'{p:.2g}')
 tex.extend([r'\begin{table}[!htbp]',r'\centering\small',r'\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{llccccc}',r'\toprule',r'Method & Control & $b/c$ & $n$ & Difference [95\% CI] & $p$ & $p_{\mathrm{Holm}}$ \\',r'\midrule'])
 for slug in generators:
  for c,cname in zip(control_slugs,['Repack','FastRelax']):
   row=r[(r.tier==1)&(r.method_a==slug)&(r.method_b==c)].iloc[0]
   effect=100*row.paired_difference;lo=100*row.cluster_bootstrap_ci_low;hi=100*row.cluster_bootstrap_ci_high
   cells=[labels[slug].replace('‡',r'$^{\ddagger}$'),cname,f'{row.a_only}/{row.b_only}',str(row.n_evaluable),f'{effect:+.1f} [{lo:+.1f}, {hi:+.1f}]',pformat(row.cluster_signflip_p),pformat(row.cluster_signflip_p_holm_controls_within_tier)]
   tex.append(' & '.join(cells)+r' \\')
 tex.extend([r'\bottomrule',r'\end{tabular}',r'\caption{\textbf{Tier 1 recovery compared with Rosetta controls.} Differences are percentage points with 95\% protein-cluster bootstrap intervals, unadjusted for multiple comparisons. Discordant pairs favor the generator ($b$) or control ($c$); $n$ counts pairs evaluable in both arms, excluding missing predictions. Two-sided protein-cluster sign-flip $p$ values are shown before and after Holm correction across the 16 prespecified comparisons of the eight generative pipelines with both Rosetta controls. A dash denotes identical paired outcomes ($b+c=0$), for which no test is performed. $\ddagger$ RFD3-generated rotamers on sampled backbones.}',r'\label{tab:recovery-controls}',r'\end{table}'])
(A/'review_statistics_control_table.tex').write_text('\n'.join(tex)+'\n')
x=s.merge(meta[['target','n_mov']],on='target',validate='many_to_one');sens=[]
for slug in slugs:
 for tier in range(1,4):
  for name,mask in [('all',np.ones(len(x),bool)),('n_mov_gt_1',x.n_mov>1)]:
   t=x[(x.method==slug)&mask];v=t[f'tier{tier}'].dropna()
   sens.append(dict(method=slug,tier=tier,subset=name,n_evaluable=len(v),n_missing=len(t)-len(v),n_recovered=int(v.sum()),recovery_fraction=v.mean()))
pd.DataFrame(sens).to_csv(A/'review_statistics_single_residue_sensitivity.csv',index=False)
plotx=x[x.method.isin(sampling_pipelines)]
bytarget=plotx.groupby('target')[['tier1','tier2','tier3']].mean().add_prefix('fraction_methods_').reset_index().merge(meta,on='target')
bytarget['n_evaluable_methods']=bytarget.target.map(plotx.groupby('target').tier1.count())
bytarget.to_csv(A/'review_statistics_target_complexity.csv',index=False)
report=['# Paired recovery statistics','',f'{len(slugs)} methods; {n_tests} planned all-pairs comparisons across three tiers.',f'Missing method-target cells (n_valid=0 or absent): {sum(v["n_missing"] for v in summ if v["tier"]==1)}.', 'Each comparison uses only pairs with valid samples for both methods; its effect denominator and protein clusters use this complete-case subset.', 'The primary family retains the eight generative pipelines versus each of two Rosetta controls (16 comparisons per tier), Holm-adjusted within tier. When included, later-added ProDy/control comparisons form a separate supporting family of six hypotheses across three tiers. ProDy is a classical backbone sampler with RFD3-generated rotamers, not a learned backbone generator. The expanded all-pairs family is a separate globally adjusted secondary analysis.', 'Identical paired outcomes have status=identical and p=NaN. They retain a planned family slot (p=1 internally for Holm adjustment), with no inferential p value reported.', 'Exactness of cluster sign flips requires joint invariance to independent method-label swaps across protein clusters; equality of mean recovery alone is insufficient.', 'For each complete-case subset, bootstrap its proteins with replacement 20,000 times, retain all pairs within each selected protein, recompute the pair-weighted effect, and take percentile 95% intervals (seed732). These intervals are pointwise, not simultaneous.',f'Secondary cluster-aware Holm-significant comparisons: {int((r.cluster_signflip_p_holm_all_pairs<.05).sum())}/{n_tests}.',f'Secondary McNemar Holm-significant comparisons: {int((r.mcnemar_exact_p_holm_all_pairs<.05).sum())}/{n_tests}.','No significant difference does not establish equivalence.']
(A/'review_statistics_report.md').write_text('\n\n'.join(report)+'\n');print('\n'.join(report));print(r.groupby(['n_evaluable','n_proteins']).size());print(pd.DataFrame(sens).query('tier==1').to_string(index=False))
if args.no_figure:raise SystemExit(0)

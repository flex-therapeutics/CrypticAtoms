from pathlib import Path
import sys,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
r=Path(__file__).parent;w=r.parent;sys.path.insert(0,str(w/'fig1'));import flexnma_plotstyle as PS
PS.apply(8)
for p in (w/'fig1/fonts').glob('*.TTF'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Verdana','font.weight':'normal','text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.edgecolor':'black','axes.linewidth':1.2,'xtick.major.width':1.1,'ytick.major.width':1.1,'pdf.fonttype':42,'svg.fonttype':'none'})
def ligand_image(path):
 im=plt.imread(path)
 mask=(im[:,:,3]>.05) if im.shape[2]==4 else np.any(im[:,:,:3]<.95,axis=2)
 ys,xs=np.nonzero(mask)
 if len(xs):
  pad=max(8,int(.025*max(xs.max()-xs.min(),ys.max()-ys.min())))
  im=im[max(0,ys.min()-pad):ys.max()+pad+1,max(0,xs.min()-pad):xs.max()+pad+1]
 return im
fig=plt.figure(figsize=(7.3,8.6));cols=['#4B7182','#A46A87'];methods=['apo2mol','molgen_ns1.0'];names=['Apo2Mol','SAGE-Flex'];targets=['KRAS','Menin','CB1','FUT8','METTL14'];rng=np.random.default_rng(15)
for row,(method,name) in enumerate(zip(methods,names)):
 top=.985-row*.31;fig.text(.013,top,'ab'[row],weight='bold',fontsize=12);fig.text(.05,top,name,fontsize=9)
 for j,t in enumerate(targets):
  ax=fig.add_axes([.012+j*.197,top-.162,.195,.13]);ax.imshow(plt.imread(r/'renders'/f'{t}_{method}.png'));ax.axis('off');ax.set_title(t,fontsize=8,pad=0)
  # Crop and fit each 2D ligand to the full target column.
  ax=fig.add_axes([.016+j*.197,top-.279,.187,.115]);ax.imshow(ligand_image(r/'renders'/f'{t}_{method}_ligand_2d.png'));ax.axis('off')
q=pd.read_csv(r/'per_target.csv');c=pd.read_csv(r/'size_motion_correlations.csv');g=pd.read_csv(r/'geometry_per_target.csv')
plots=[('bond_issue_percent','Ligands with bond /\nconnectivity issues (%)',g,None),('intraligand_clash_percent','Within-ligand clashes\n(% ligand atoms)',g,None),('protein_ligand_clash_percent','Protein–ligand clashes\n(% ligand atoms)',g,None),('tanimoto','Reference Tanimoto similarity',q,None),('rg_ratio','Ligand radius / reference',q,1),('spearman_rg_Dsa','Ligand size vs. pocket motion (ρ)',c,0)]
for k,(metric,label,data,baseline) in enumerate(plots):
 x=.077+(k%3)*.322;y=.215 if k<3 else .035;ax=fig.add_axes([x,y,.235,.105]);ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3,labelsize=7);fig.text(x-.047,y+.139,'cdefgh'[k],fontsize=12,weight='bold');ax.set_title(label,fontsize=7.3,pad=9)
 if baseline is not None:ax.axhline(baseline,c='#888888',ls=':',lw=.8,zorder=0)
 for i,m in enumerate(methods):
  vals=data.loc[data.method==m,metric].dropna().values;ax.scatter(i+rng.uniform(-.13,.13,len(vals)),vals,s=12,c=cols[i],alpha=.4,edgecolors='none',zorder=1)
  ax.boxplot([vals],positions=[i],widths=.42,patch_artist=True,showfliers=False,manage_ticks=False,boxprops={'facecolor':'none','edgecolor':cols[i],'linewidth':1.3},medianprops={'color':cols[i],'linewidth':1.8},whiskerprops={'color':cols[i]},capprops={'color':cols[i]},zorder=3)
 ax.set_xticks([0,1],names,fontsize=7);ax.set_xlim(-.6,1.6)
 if metric=='tanimoto':ax.set_ylim(0,.2)
 if metric.endswith('percent'):ax.set_ylim(bottom=-1)
 if metric=='spearman_rg_Dsa':ax.set_ylim(-.65,.65)
for ext in ['pdf','svg','png']:fig.savefig(r/f'figureA11.{ext}',dpi=300,transparent=True)

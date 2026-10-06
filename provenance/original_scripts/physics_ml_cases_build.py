from pathlib import Path
import sys,json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
r=Path(__file__).parent;w=r.parent;sys.path.insert(0,str(w/'fig1'));import flexnma_plotstyle as PS
PS.apply(8)
for p in (w/'fig1/fonts').glob('*.TTF'):font_manager.fontManager.addfont(str(p))
plt.rcParams.update({'font.family':'Verdana','font.weight':'normal','text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.edgecolor':'black','axes.linewidth':1.2,'xtick.major.width':1.1,'ytick.major.width':1.1,'mathtext.fontset':'custom','mathtext.rm':'Verdana','mathtext.it':'Verdana:italic','mathtext.bf':'Verdana:bold','pdf.fonttype':42,'svg.fonttype':'none'})
d=pd.read_csv(r/'per_sample.csv');refs=pd.read_csv(r/'references.csv').set_index('target');sel=pd.read_csv(r/'selected_examples.csv');cfg={x[0]:(x[1],x[3]) for x in json.loads((w/'controls/methods_with_controls.json').read_text())};methods=['rosetta_repack','rosetta_fastrelax','prody_rfd3','bioemu_refine','boltz2','esmfold2'];marks=['s','^','v','o','P','D'];fig=plt.figure(figsize=(7.3,7.7))
for x,title in zip([.112,.332,.552,.844],['Apo–holo overlay','Closest physics sample','Recovered ML sample','Backbone and side-chain error']):fig.text(x,.965,title,ha='center',fontsize=7.1)
for i,t in enumerate(['cAbl','HIV_RT__3V81','GLTP']):
 top=.912-i*.282;fig.text(.008,top+.018,'abc'[i],fontsize=12,weight='bold');fig.text(.042,top+.018,{'cAbl':'c-Abl','HIV_RT__3V81':'HIV RT (3V81)'}.get(t,t),fontsize=9)
 for j,group in enumerate(['apo','physics','learned']):
  ax=fig.add_axes([.001+j*.22,top-.191,.22,.195]);ax.imshow(plt.imread(r/'renders'/f'{t}_{group}.png'));ax.axis('off')
  if group=='apo':line=rf"$D_{{\mathrm{{ah}}}}$ = {refs.loc[t,'D']:.2f} Å";name='Apo (gray), holo (teal)'
  else:
   q=sel[sel.target.eq(t)&sel.group.eq(group)].iloc[0];name=cfg[q.method][0];line=rf'$D_{{\mathrm{{sa}}}}$ {q.Dsa:.2f} · $D_{{\mathrm{{sh}}}}$ {q.Dsh:.2f} Å'
  fig.text(.111+j*.22,top-.203,name,fontsize=7,ha='center');fig.text(.111+j*.22,top-.222,line,fontsize=6.7,ha='center')
 ax=fig.add_axes([.727,top-.19,.255,.185]);ax.spines[['top','right']].set_visible(False);q=d[d.target.eq(t)]
 for m,mk in zip(methods,marks):
  v=q[q.method.eq(m)];ax.scatter(v.bb_rmsd,v.sc_rmsd,c=cfg[m][1],marker=mk,s=13,alpha=.45,edgecolors='none',zorder=2)
 ax.scatter([refs.loc[t,'bb_rmsd']],[refs.loc[t,'sc_rmsd']],marker='*',s=55,c='black',zorder=5)
 for row in sel[sel.target.eq(t)].itertuples():ax.scatter([row.bb_rmsd],[row.sc_rmsd],s=50,facecolors='none',edgecolors='black',lw=.8,zorder=6)
 ax.set_xlim(left=0);ax.set_ylim(bottom=0);ax.tick_params(labelsize=6.5,length=3,pad=2);ax.set_xlabel('Backbone RMSD (Å)',fontsize=7,labelpad=2);ax.set_ylabel('Side-chain RMSD (Å)',fontsize=7,labelpad=2)
handles=[Line2D([],[],marker=mk,ls='none',color=cfg[m][1],label=cfg[m][0],markersize=4) for m,mk in zip(methods,marks)]+[Line2D([],[],marker='*',ls='none',color='black',label='Apo reference',markersize=6)]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.51,.003),ncol=4,frameon=False,fontsize=6.8,columnspacing=1.2,handletextpad=.4)
for ext in ['pdf','svg','png']:fig.savefig(r/f'figureA12.{ext}',dpi=300,transparent=True)

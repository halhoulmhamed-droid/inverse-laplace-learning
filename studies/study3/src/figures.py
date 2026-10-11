"""Figures from archived study-3 losses and selections, no experiment rerun."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import NullLocator
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'figures';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader((BASE/'results/analysis/summary.csv').open()))
report=json.loads((BASE/'results/analysis/report.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
manifest=[]
def save(name,fig,**arrays):
 np.savez_compressed(OUT/(name+'.npz'),**arrays);fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
 manifest.append({'name':name,'data':f'figures/{name}.npz','figure':f'figures/{name}.pdf','script':'src/figures.py','configuration':'results/final/configuration.json','analysis':'src/analysis.py'})
def values(m,scenario='primary',N=64,key='risk'):
 return np.array([float(next(x for x in rows if x['method']==m and x['scenario']==scenario and x['N']==str(N) and x['game']==str(g))[key]) for g in range(5)])
# Primary and protection: upper bounds only, same physical-object bootstrap.
fig,axes=plt.subplots(1,2,figsize=(9.6,3.4),layout='constrained');p=report['primary']
data={k:np.array([x[k] for x in p]) for k in ['primary_ratio','primary_upper','protection_ratio','protection_upper']}
for ax,prefix,limit,title in zip(axes,['primary','protection'],[.97,1.03],['Gain utile : capteurs complets','Protection : 32 capteurs contigus']):
 for i in range(5):ax.plot([data[prefix+'_ratio'][i],data[prefix+'_upper'][i]],[i]*2,color='#17747e',lw=2)
 ax.scatter(data[prefix+'_ratio'],range(5),color='#193049',s=30,label='Ratio observé',zorder=3)
 ax.scatter(data[prefix+'_upper'],range(5),color='#17747e',marker='|',s=130,label='Borne supérieure 99,5 %',zorder=3)
 ax.axvline(limit,color='#b25b25',ls='--',label=f'Seuil {limit:.2f}');ax.axvline(1,color='gray',lw=.8)
 ax.set(xlim=(.945,1.04),yticks=range(5),xlabel='Risque S / M0',ylabel='Jeu d’entraînement',title=title);ax.grid(alpha=.2)
axes[0].legend(fontsize=7,loc='lower right')
save('primary',fig,**data,thresholds=np.array([.97,1.03]))
# Training sizes are descriptive, including N32 and256, not selected post hoc.
fig,axes=plt.subplots(1,2,figsize=(9.4,3.3),layout='constrained');sizes=np.array([32,64,256]);arrays={}
for ax,scenario,title in zip(axes,['primary','partial'],['Capteurs complets','32 capteurs contigus']):
 ratios=np.array([values('S',scenario,n)/values('M0',scenario,n) for n in sizes]);arrays[scenario]=ratios
 for g in range(5):ax.plot(sizes,ratios[:,g],'-o',alpha=.55,lw=.9,ms=4)
 ax.plot(sizes,ratios.mean(1),'k-',lw=2,label='Moyenne de cinq ratios');ax.axhline(1,color='gray',ls='--');ax.set(xscale='log',xticks=sizes,xticklabels=sizes,ylabel='Risque S / M0',xlabel="Objets d'apprentissage",title=title);ax.grid(alpha=.2)
for ax in axes: ax.xaxis.set_minor_locator(NullLocator())
axes[0].legend(fontsize=8)
save('training',fig,sizes=sizes,**arrays)
scenarios=['primary','partial','even32','quarter','noise_low','noise_high','rough_full','rough_partial','height_full','height_partial']
labels=['Complet, nominal','32 contigus, nominal','32 espacés','16 contigus','Complet, bruit 0,003','Complet, bruit 0,03','Complet, rugueux','32 contigus, rugueux','Complet, hauteur +10 %','32 contigus, hauteur +10 %']
pairs=[('S','M0'),('S','E'),('S','H1'),('Diag','M0')]
ratios=np.array([[values(a,s).mean()/values(b,s).mean() for a,b in pairs] for s in scenarios])
fig,ax=plt.subplots(figsize=(8.7,4.8),layout='constrained');im=ax.imshow(np.log2(ratios),cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-.7,vcenter=0,vmax=.15),aspect='auto')
ax.set(yticks=range(10),yticklabels=labels,xticks=range(4),xticklabels=['S / M0','S / E','S / H1','Diag / M0'],title='N=64 : ratios de risques moyens sur cinq entraînements')
for i in range(10):
 for j in range(4):ax.text(j,i,f'{ratios[i,j]:.3f}',ha='center',va='center',color='white' if ratios[i,j]<.75 else 'black')
fig.colorbar(im,ax=ax,label='log₂(ratio) ; < 0 : perte plus faible',shrink=.8)
save('robustness',fig,ratios=ratios,scenarios=np.array(scenarios),pairs=np.array(pairs))
# Display actual gamma values, labels discrete for irregular grid.
keys=[('full',.01),('half',.01),('even32',.01),('quarter',.01),('full',.003),('full',.03)]
grid=np.array([0,.01,.03,.1,.3,.6,1]);gammas=np.empty((15,6));alphas=np.empty_like(gammas)
for i,(N,g) in enumerate((n,g) for n in sizes for g in range(5)):
 for j,(m,noise) in enumerate(keys):
  selected=json.loads((BASE/f'results/final/n{N}_g{g}_{m}_{noise}_selection.json').read_text())['selected']['S'];gammas[i,j]=selected['gamma'];alphas[i,j]=selected['parameter']
indices=np.searchsorted(grid,gammas)
fig,ax=plt.subplots(figsize=(8.7,5.1),layout='constrained');im=ax.imshow(indices,cmap='YlGnBu',vmin=0,vmax=6,aspect='auto')
ax.set(xticks=range(6),xticklabels=['Complet','32 contigus','32 espacés','16 contigus','Bruit 0,003','Bruit 0,03'],yticks=range(15),yticklabels=[f'N={n} · jeu {g}' for n in sizes for g in range(5)],title='Poids γ sélectionné uniquement sur validation nominale')
for i in range(15):
 for j in range(6):ax.text(j,i,f'{gammas[i,j]:g}',ha='center',va='center',color='white' if indices[i,j]>=4 else 'black',fontsize=8)
for y in [4.5,9.5]:ax.axhline(y,color='white',lw=2)
save('gamma',fig,gammas=gammas,alphas=alphas,gamma_grid=grid,sizes=sizes,keys=np.array([f'{a}_{b}' for a,b in keys]))
# Matched alpha compared with independently tuned baseline; no causal overclaim.
methods=['S','M0','M0_at_S','H1m','H1m_at_S','E','H1','Diag'];means={}
fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for ax,scenario,title in zip(axes,['primary','partial'],['Complet, N=64','32 contigus, N=64']):
 v=np.array([values(m,scenario).mean() for m in methods]);means[scenario]=v
 ax.bar(np.arange(8),v,color=['#17747e','#55778d','#d79d49','#4f9e8f','#d9b979','#888888','#a67b90','#718b50'])
 ax.set(xticks=range(8),xticklabels=['S','M0','M0 | αS','H1m','H1m | αS','E','H1','Diag'],ylabel='Risque quadratique L²',title=title);ax.tick_params(axis='x',labelrotation=35);ax.grid(axis='y',alpha=.2)
save('matched',fig,methods=np.array(methods),**means)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
# Reproducible tables read by both LaTeX documents.
docs=BASE/'docs';docs.mkdir(exist_ok=True)
lines=[r'\begin{tabular}{rrrrr}',r'\toprule Jeu & Ratio complet & Borne sup. & Ratio contigu & Borne sup.\\\midrule']
for x in p:lines.append(str(x['game'])+' & '+' & '.join(f"{x[k]:.4f}" for k in ['primary_ratio','primary_upper','protection_ratio','protection_upper'])+r'\\')
lines.extend([r'\bottomrule',r'\end{tabular}']);(docs/'primary_table.tex').write_text('\n'.join(lines)+'\n')
lines=[r'\begin{tabular}{lrr}',r'\toprule Méthode & Complet & 32 contigus\\\midrule']
for m in ['S','M0','Diag','E','H1','H1m']:lines.append(m+f" & {values(m).mean():.7f} & {values(m,'partial').mean():.6f}"+r'\\')
lines.extend([r'\bottomrule',r'\end{tabular}']);(docs/'method_table.tex').write_text('\n'.join(lines)+'\n')
print(json.dumps({'figures':len(manifest),'gamma_zero':int(np.sum(gammas==0)),'gamma_one':int(np.sum(gammas==1)),'total_fits':gammas.size}))

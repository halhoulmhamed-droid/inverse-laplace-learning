"""Render frozen summaries; every figure gets its exact numerical NPZ."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'figures';OUT.mkdir(exist_ok=True)
rows=list(csv.DictReader((BASE/'results/analysis/summary.csv').open()))
report=json.loads((BASE/'results/analysis/report.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
manifest=[]
def save(name,fig,**arrays):
 np.savez_compressed(OUT/(name+'.npz'),**arrays)
 fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
 manifest.append(dict(name=name,data=f'figures/{name}.npz',figure=f'figures/{name}.pdf',script='src/figures.py',analysis='src/analysis.py',configuration='results/final/configuration.json'))
def values(method,scenario='primary',N=256,metric='risk'):
 return np.array([float(next(x for x in rows if x['method']==method and x['scenario']==scenario and int(x['N'])==N and int(x['game'])==g)[metric]) for g in range(5)])
# Conditional uncertainty, no population training CI.
f,axes=plt.subplots(1,2,figsize=(9.4,3.4),layout='constrained')
ps=report['primary'];data={k:np.array([p[k] for p in ps]) for k in ['risk_ratio','risk_upper','tail_ratio','tail_upper']}
for ax,point,upper,threshold,title in zip(axes,['risk_ratio','tail_ratio'],['risk_upper','tail_upper'],[.95,1.05],['Risque moyen F / D','95e percentile F / D']):
 for i in range(5):ax.plot([data[point][i],data[upper][i]],[i+1,i+1],color='#17747e',lw=2)
 ax.scatter(data[point],np.arange(1,6),color='#193049',label='Ratio observé',zorder=4)
 ax.scatter(data[upper],np.arange(1,6),marker='|',s=130,color='#17747e',label='Borne supérieure 99,5 %')
 ax.axvline(threshold,color='#b65d27',ls='--',label=f'Seuil {threshold:.2f}');ax.set(xlim=(.2,1.12),yticks=range(1,6),ylabel='Entraînement',title=title,xlabel='Ratio de pertes quadratiques');ax.grid(alpha=.2)
axes[0].legend(fontsize=7,loc='upper right')
save('primary',f,**data,training=np.arange(5),threshold=np.array([.95,1.05]))
# Learning size, separate across-training spread.
f,axes=plt.subplots(1,2,figsize=(9.4,3.4),layout='constrained');sizes=np.array([64,256,1024]);d={}
for ax,a,b in zip(axes,['F','L'],['D','H1']):
 rat=np.array([values(a,N=n)/values(b,N=n) for n in sizes]);d[a+'_over_'+b]=rat
 for g in range(5):ax.plot(sizes,rat[:,g],'-o',alpha=.45,lw=.9,ms=4)
 ax.plot(sizes,rat.mean(1),'k-',lw=2,label='Moyenne des 5 ratios');ax.axhline(1,color='gray',ls='--');ax.set(xscale='log',xticks=sizes,xticklabels=sizes,xlabel="Objets d'entraînement",ylabel=f'Risque {a} / {b}',title=f'{a} face à {b} — scénario principal');ax.grid(alpha=.2)
axes[1].legend(fontsize=8)
save('training',f,sizes=sizes,**d)
# Robustness: ratio of mean risks, descriptive.
scenarios=['primary','full','even32','quarter','even16','noise_low','noise_high','rough','height5','height10']
labels=['32 contigus, nominal','64 complets','32 espacés','16 contigus','16 espacés','Bruit 0,003','Bruit 0,03','Rugosité p=1,5','Hauteur +5 %','Hauteur +10 %']
pairs=[('F','D'),('L','H1'),('O','F')]
ratios=np.array([[values(a,s).mean()/values(b,s).mean() for a,b in pairs] for s in scenarios])
f,ax=plt.subplots(figsize=(7.8,4.6),layout='constrained');im=ax.imshow(np.log2(ratios),cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=1),aspect='auto')
ax.set(yticks=range(10),yticklabels=labels,xticks=range(3),xticklabels=['F / D','L / H1','O / F'],title='N=256 : ratios de risques moyens sur cinq entraînements')
for i in range(10):
 for j in range(3):ax.text(j,i,f'{ratios[i,j]:.3f}',ha='center',va='center',color='black' if ratios[i,j]>.45 else 'white',fontsize=10)
f.colorbar(im,ax=ax,label='log₂(ratio) ; < 0 : perte plus faible',shrink=.8)
save('robustness',f,ratios=ratios,scenarios=np.array(scenarios),pairs=np.array(pairs))
# Shared-ridge decomposition: full/projected errors identify complementary prediction.
methods=['D_at_F','F_at_F','O_at_F','L_at_L','LP_at_L','H1']
projected=np.array([values(m,metric='projected').mean() for m in methods]);comp=np.array([values(m,metric='complementary').mean() for m in methods])
f,ax=plt.subplots(figsize=(9.2,3.8),layout='constrained');x=np.arange(6)
ax.bar(x,projected,color='#17747e',label='Erreur dans P');ax.bar(x,comp,bottom=projected,color='#d59b44',label='Complément + modes 64–127')
ax.set(xticks=x,xticklabels=['D | λF','F | λF','O | λF','L | λL','LP | λL','H1 réglé'],ylabel='Risque quadratique L²',title='Décomposition à régularisation partagée — N=256, 32 contigus');ax.legend(fontsize=9);ax.grid(axis='y',alpha=.2)
save('projection',f,methods=np.array(methods),projected=projected,complementary=comp)
# Analytic population covariance in each stored numerical SVD basis.
f,axes=plt.subplots(1,3,figsize=(10,3.4),layout='constrained');covs=[];fractions=[]
w=(1+np.arange(128,dtype=float)**2)**(-2.5);w/=w.sum()
for ax,mask,label in zip(axes,['full','half','quarter'],['64 complets','32 contigus','16 contigus']):
 V=np.load(BASE/f'results/final/n256_g0_{mask}_0.01_models.npz')['V'];C=V.T@(w[:64,None]*V)
 covs.append(C);fraction=np.linalg.norm(C-np.diag(np.diag(C)))/np.linalg.norm(C);fractions.append(fraction)
 im=ax.imshow(C/max(np.diag(C)),vmin=-1,vmax=1,cmap='RdBu_r');ax.set(title=f'{label}\nhors diagonale : {fraction:.1%}',xlabel='Direction singulière',ylabel='Direction singulière',xticks=[0,5,11],xticklabels=[1,6,12],yticks=[0,5,11],yticklabels=[1,6,12])
f.colorbar(im,ax=axes,label='Covariance / variance maximale',shrink=.7)
save('covariance',f,covariance=np.array(covs),off_diagonal_fraction=np.array(fractions),cosine_variances=w)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'figures':len(manifest),'population_off_diagonal_fractions':fractions}))

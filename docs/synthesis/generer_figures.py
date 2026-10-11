"""Mise en graphique de valeurs publiees dans les rapports E2 et E3.

Aucune simulation, estimation ou inference nouvelle. Les CSV sont des
transcriptions des tableaux arrondis ; ils ne remplacent pas les pertes brutes.
"""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'axes.edgecolor':'#9ba8b1','text.color':'#173044',
                     'axes.labelcolor':'#173044','xtick.color':'#425563',
                     'ytick.color':'#425563','pdf.fonttype':42})
with (ROOT/'data/etude2_risques.csv').open() as f:
    rows=list(csv.DictReader(f))
fig, ax = plt.subplots(figsize=(7.1,2.35), layout='constrained')
values=[float(r['risk']) for r in rows]
labels=['D · diagonal (rang 12)','F · complet (rang 12)',
        'H¹ · référence classique','L · LMMSE complet']
ax.barh(range(4),values,color=['#c4cdd4','#087e8b','#536978','#ad8444'],height=.57)
ax.set_yticks(range(4),labels); ax.invert_yaxis()
ax.set_xlim(0,.186); ax.set_xticks([0,.04,.08,.12,.16])
ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:.2f}'.replace('.',',')))
ax.set_xlabel('Risque quadratique L² moyen (plus petit = meilleur)')
ax.grid(axis='x',alpha=.15); ax.set_axisbelow(True)
for i,v in enumerate(values): ax.text(v+.0025,i,f'{v:.6f}'.replace('.',','),va='center',fontsize=9)
fig.savefig(ROOT/'figures/etude2_risques.pdf',metadata={'Title':'Etude 2 - risques rapportes','Author':'Halhoul Mhamed','CreationDate':None});plt.close(fig)

with (ROOT/'data/etude3_criteres.csv').open() as f: rows=list(csv.DictReader(f))
fig,axes=plt.subplots(1,2,figsize=(7.1,2.9),sharey=True,layout='constrained')
for ax,key,title,threshold in zip(axes,['full','contiguous'],['Gain utile · 64 capteurs','Protection · 32 contigus'],[.97,1.03]):
    for i,r in enumerate(rows):
        x=float(r[key+'_ratio']);u=float(r[key+'_upper'])
        ax.plot([x,u],[i,i],color='#9caeb9',lw=1.7,zorder=1)
        ax.scatter(x,i,color='#087e8b',s=32,zorder=3)
        ax.scatter(u,i,color='#ad8444',marker='|',s=100,zorder=4)
    ax.axvline(threshold,color='#a54d44',ls='--',lw=1.1)
    ax.axvline(1,color='#aab3ba',lw=.7)
    ax.set_xlim(.951,1.038);ax.set_xticks([.96,.98,1,1.02])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:.2f}'.replace('.',',')))
    ax.set_title(title,fontsize=11,pad=12);ax.set_xlabel('Ratio de risques S / M0')
    ax.set_yticks(range(5),[f'Jeu {i}' for i in range(5)])
    ax.grid(axis='x',alpha=.12)
    ax.text(threshold,.98,f'Seuil {threshold:.2f}'.replace('.',','),transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=8,color='#a54d44')
axes[0].invert_yaxis()
fig.legend(handles=[plt.Line2D([],[],marker='o',color='#087e8b',ls='',label='Ratio observé'),plt.Line2D([],[],marker='|',color='#ad8444',ls='',markersize=10,label='Borne supérieure 99,5 %')],loc='outside lower center',ncol=2,frameon=False,fontsize=9)
fig.savefig(ROOT/'figures/etude3_criteres.pdf',metadata={'Title':'Etude 3 - criteres rapportes','Author':'Halhoul Mhamed','CreationDate':None});plt.close(fig)
print('Deux figures produites à partir des seuls CSV transcrits.')

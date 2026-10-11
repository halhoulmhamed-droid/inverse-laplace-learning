"""Projection layout of the preselected P1 worst object; no new experiment."""
from pathlib import Path
import hashlib,json,resource
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
source=root/'figures/P1_pire.npz'
with np.load(source) as z: data=dict(z)
x=data['x'];truth=data['truth_coefficients'];c=data['reconstructed_coefficients']
phi=np.sqrt(2)*np.cos(np.pi*x[:,None]*np.arange(len(truth)));phi[:,0]=1
names=['Tikhonov L2','Tikhonov H1','TSVD','Wiener','Wiener plafonné','Appris signé','Signé plafonné','Appris positif','Appris borné','Borné plafonné','LMMSE']
plt.rcParams.update({'font.size':11,'pdf.fonttype':42})
fig,axes=plt.subplots(1,3,figsize=(12,3.5),sharex=True,sharey=True)
for ax,indices in zip(axes,[range(4),range(4,8),range(8,11)]):
 ax.plot(x,phi@truth,color='black',lw=2.5,label='Vérité')
 for j in indices: ax.plot(x,phi[:,:c.shape[1]]@c[j],lw=1.5,ls=['-','--',':','-.'][j%4],label=names[j])
 ax.set(xlabel='x',ylabel='q(x)');ax.legend(fontsize=9,loc='upper right');ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for ext in ['pdf','png']:fig.savefig(root/('figures/P1_pire_projection.'+ext),dpi=140,bbox_inches='tight')
record={'source':'figures/P1_pire.npz','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'script':'scripts/presentation_figures.py','object':197,'repetition':0,'methods':'all eleven, same inputs','cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime}
(root/'figures/P1_pire_projection.json').write_text(json.dumps(record,indent=2)+'\n')

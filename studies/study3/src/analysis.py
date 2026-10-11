"""Predefined risk/protection inference, paired jointly by physical object."""
import csv,json,hashlib,time
from pathlib import Path
import numpy as np
from studies.study3.src.data import rng
from studies.study3.src.model import METHODS,FIXED
BASE=Path(__file__).resolve().parents[1]
def bootstrap(primary,partial,game,draws=4000):
 # Both inputs: object x repeat x two methods [S,M0].
 a=np.mean(primary,axis=1);b=np.mean(partial,axis=1)
 if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=2:raise ValueError('paired dimensions')
 if not all(np.isfinite(x).all() for x in [a,b]) or a[:,1].mean()<=0 or b[:,1].mean()<=0:
  return {'valid':False,'reason':'invalid input or denominator'},np.full((draws,2),np.nan)
 random=rng(40,game,0);vals=[]
 for start in range(0,draws,200):
  idx=random.integers(0,len(a),size=(min(200,draws-start),len(a)))
  means_a=a[idx].mean(1);means_b=b[idx].mean(1)
  vals.append(np.column_stack([means_a[:,0]/means_a[:,1],means_b[:,0]/means_b[:,1]]))
 vals=np.concatenate(vals)
 if not np.isfinite(vals).all():return {'valid':False,'reason':'nonfinite bootstrap'},vals
 upper=np.quantile(vals,.995,axis=0)
 r=dict(valid=True,primary_ratio=float(a[:,0].mean()/a[:,1].mean()),protection_ratio=float(b[:,0].mean()/b[:,1].mean()),primary_upper=float(upper[0]),protection_upper=float(upper[1]),draws=draws,upper_quantile=.995)
 r['criteria_met']=bool(upper[0]<=.97 and upper[1]<=1.03)
 return r,vals

def main():
 t=time.process_time();f=BASE/'results/final';out=BASE/'results/analysis'
 if out.exists():raise RuntimeError('Existing analysis preserved')
 out.mkdir();cfg=json.loads((f/'configuration.json').read_text());rows=[];decisions=[];boots=[];failures=0
 for N in cfg['sizes']:
  for g in cfg['games']:
   conditions={}
   for name,*_ in cfg['scenarios']:
    z=np.load(f/f'n{N}_g{g}_{name}_loss.npz');loss=z['loss'];status=z['status'];failures+=int(np.count_nonzero(status));conditions[name]=loss
    for j,m in enumerate(METHODS+FIXED):
     valid=not np.any(status[:,:,j]) and np.isfinite(loss[:,:,j]).all()
     row=dict(N=N,game=g,scenario=name,method=m,valid=bool(valid),failures=int(np.count_nonzero(status[:,:,j])))
     for k,n in enumerate(['risk','high_mode','residual']):row[n]=float(loss[:,:,j,k].mean()) if valid else None
     row['q95']=float(np.quantile(loss[:,:,j,0].mean(1),.95)) if valid else None;rows.append(row)
   if N==64:
    result,v=bootstrap(conditions['primary'][:,:,:2,0],conditions['partial'][:,:,:2,0],g)
    decisions.append(dict(game=g,**result));boots.append(v)
 with (out/'summary.csv').open('w') as s:
  w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 np.savez_compressed(out/'bootstrap.npz',ratios=np.stack(boots))
 result=dict(primary=decisions,conclusion='ESTABLISHED_ON_FROZEN_CONDITIONS' if all(x.get('valid') and x.get('criteria_met') for x in decisions) else 'NOT_ESTABLISHED',
  rows=len(rows),failures=failures,cpu_seconds=time.process_time()-t,
  raw_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in f.glob('*_loss.npz')},
  inference='Paired object percentile bootstrap, same indices for full/partial, 99.5% upper bounds for 5 training runs x 2 inequalities; approximate conditional coverage, no population training confidence claim.')
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='raw_sha256'},indent=2))
if __name__=='__main__':main()

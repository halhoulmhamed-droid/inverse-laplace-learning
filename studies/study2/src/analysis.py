"""Frozen primary inference; no tuning or data exclusion."""
import csv,json,hashlib,time
from pathlib import Path
import numpy as np
from studies.study2.src.model import METHODS,FIXED,rng
BASE=Path(__file__).resolve().parents[1]

def paired_bootstrap(loss_a,loss_b,game,draws=4000):
    # Inputs object x repeats. Never resample noise repetitions separately.
    a,b=np.mean(loss_a,axis=1),np.mean(loss_b,axis=1)
    if a.shape!=b.shape or a.ndim!=1: raise ValueError('paired object dimensions')
    if not np.isfinite(a).all() or not np.isfinite(b).all() or b.mean()<=0 or np.quantile(b,.95)<=0:
        return dict(valid=False,reason='invalid or zero denominator'),np.full((draws,2),np.nan)
    random=rng(40,game,0);values=[]
    for start in range(0,draws,200):
        idx=random.integers(0,len(a),size=(min(200,draws-start),len(a)))
        aa,bb=a[idx],b[idx]
        values.append(np.column_stack([aa.mean(1)/bb.mean(1),np.quantile(aa,.95,axis=1)/np.quantile(bb,.95,axis=1)]))
    boot=np.concatenate(values)
    if not np.isfinite(boot).all(): return dict(valid=False,reason='invalid bootstrap'),boot
    upper=np.quantile(boot,.995,axis=0)
    result=dict(valid=True,risk_ratio=float(a.mean()/b.mean()),tail_ratio=float(np.quantile(a,.95)/np.quantile(b,.95)),
                risk_upper=float(upper[0]),tail_upper=float(upper[1]),draws=draws,quantile=.995)
    result['criteria_met']=bool(upper[0]<=.95 and upper[1]<=1.05)
    return result,boot

def main():
    start=time.process_time();folder=BASE/'results/final';cfg=json.loads((folder/'configuration.json').read_text())
    out=BASE/'results/analysis'
    if out.exists(): raise RuntimeError('Preserve previous analysis')
    out.mkdir()
    records=[];primary=[];boots=[];invalid=0
    for N in cfg['sizes']:
        for g in cfg['games']:
            for s,*_ in cfg['scenarios']:
                with np.load(folder/f'n{N}_g{g}_{s}_loss.npz') as z:
                    loss=z['loss'];status=z['status']
                invalid+=int(np.count_nonzero(status))
                for j,m in enumerate(METHODS+FIXED):
                    entry=dict(N=N,game=g,scenario=s,method=m,failures=int(np.count_nonzero(status[:,:,j])))
                    valid=entry['failures']==0 and np.isfinite(loss[:,:,j]).all()
                    entry['valid']=valid
                    for k,name in enumerate(['risk','projected','complementary']):
                        entry[name]=float(loss[:,:,j,k].mean()) if valid else None
                    entry['q95']=float(np.quantile(loss[:,:,j,0].mean(1),.95)) if valid else None
                    records.append(entry)
                if N==256 and s=='primary':
                    result,boot=paired_bootstrap(loss[:,:,1,0],loss[:,:,0,0],g)
                    primary.append(dict(game=g,**result));boots.append(boot)
    with (out/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    np.savez_compressed(out/'bootstrap.npz',ratios=np.stack(boots))
    verdict=all(p.get('valid') and p.get('criteria_met') for p in primary)
    answer=dict(primary=primary,conclusion='ESTABLISHED_ON_FROZEN_CONDITIONS' if verdict else 'NOT_ESTABLISHED',
        total_failure_entries=invalid,rows=len(records),cpu_seconds=time.process_time()-start,
        uncertainty='Conditional paired object bootstrap. Approximate family-wise 95% upper bounds across 5 trainings x 2 ratios. No population training confidence interval.',
        raw_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.glob('*_loss.npz'))})
    (out/'report.json').write_text(json.dumps(answer,indent=2)+'\n')
    print(json.dumps({k:v for k,v in answer.items() if k!='raw_sha256'},indent=2))
if __name__=='__main__': main()

"""Frozen pilot/final campaign, invoked by bounded.py from repository root."""
import argparse,hashlib,json,time,resource,sys
from pathlib import Path
import numpy as np
from laplace_learning.direct import fd_matrix
from studies.study2.src.model import *

BASE=Path(__file__).resolve().parents[1]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj): p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def meta(c): return {k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in c.items() if k!='R'}

def run(pilot=False):
    start=time.process_time();wall=time.monotonic()
    frozen=json.loads((BASE/'FREEZE.json').read_text())
    for p,h in frozen['sha256'].items(): assert sha(BASE/p)==h
    test=json.loads((BASE/'results/tests.json').read_text());assert test['success']
    for p,h in test['source_hashes'].items(): assert sha(BASE/p)==h
    out=BASE/'results'/('pilot' if pilot else 'final')
    if out.exists(): raise RuntimeError('Existing campaign preserved: refusing overwrite')
    if not pilot:
        prior=json.loads((BASE/'results/pilot/report.json').read_text())
        assert prior['success'] and prior['feasibility_pass']
    out.mkdir(parents=True)
    ns=[64] if pilot else [64,256,1024];games=range(1 if pilot else 5)
    nv,nt,reps=(24,24,2) if pilot else (256,400,4)
    scenarios=[SCENARIOS[i] for i in [0,1,3,9]] if pilot else SCENARIOS
    offset=100 if pilot else 0
    cfg=dict(pilot=pilot,root_seed=ROOT,sizes=ns,games=list(games),validation=nv,test=nt,repetitions=reps,
             scenarios=scenarios,methods=METHODS,fixed_methods=FIXED,loss_axes=['total','projected','complementary'],
             rank=12,protocol_sha256=sha(BASE/'PROTOCOL.md'),freeze_commit='20a663d',
             sources={str(p.relative_to(BASE)):sha(p) for p in (BASE/'src').glob('*.py')},
             reused_direct_sha256=sha(BASE.parents[1]/'laplace_learning/direct.py'),numpy=np.__version__,python=sys.version)
    dump(out/'configuration.json',cfg)
    forward=fd_matrix(64);failures=[];diagnostics=[];max_identity=0.;reconstructions=0
    truth={2.5:profiles(30+offset,0,nt),1.5:profiles(31+offset,0,nt,p=1.5)}
    eps={p:noises(part+offset,0,nt,reps) for p,part in [(2.5,30),(1.5,31)]}
    np.savez_compressed(out/'truth.npz',nominal=truth[2.5],rough=truth[1.5])
    for game in games:
        train=profiles(10+offset,game,max(ns));tn=noises(10+offset,game,max(ns),1)
        val=profiles(20+offset,game,nv);vn=noises(20+offset,game,nv,2)
        for N in ns:
            fitted={}
            for scenario,m,sigma,p,height in scenarios:
                idx=mask(m);key=(m,sigma)
                if key not in fitted:
                    A=forward[idx]/8
                    y=observations(train[:N],tn[:N],idx,sigma)[:,0]
                    family=covariance_family(A,train[:N,:64],y)
                    yv=observations(val,vn,idx,sigma).reshape(-1,len(idx))
                    selected,fixed,scores=select(family,yv,np.repeat(val,2,axis=0))
                    fitted[key]=(family,selected,fixed)
                    label=f'n{N}_g{game}_{m}_{sigma}'
                    dump(out/(label+'_selection.json'),dict(selected={k:meta(v) for k,v in selected.items()},candidates=[meta(a) for a in scores]))
                    np.savez_compressed(out/(label+'_models.npz'),R=np.stack([selected[k]['R'] for k in METHODS]),
                        fixed_R=np.stack([a['R'] for a in fixed]),mc=family['mc'],my=family['my'],
                        U=family['u'],V=family['v'],Cbz=family['Cbz'],Czz=family['Czz'],singular=family['s'])
                    ident=max(np.max(abs(family['P']@fixed[2]['R']-fixed[1]['R'])),
                              np.max(abs(family['P']@fixed[3]['R']-fixed[4]['R'])))
                    max_identity=max(max_identity,float(ident))
                    if ident>1e-9: failures.append(label+': projection identity')
                family,selected,fixed=fitted[key]
                yt=observations(truth[p],eps[p],idx,sigma,height).reshape(-1,len(idx))
                ct=np.repeat(truth[p],reps,axis=0)
                all_candidates=[selected[k] for k in METHODS]+fixed
                loss=np.empty((nt,reps,len(all_candidates),3));status=np.zeros((nt,reps,len(all_candidates)),dtype=np.uint8)
                for j,a in enumerate(all_candidates):
                    if not a['valid']:
                        loss[:,:,j]=np.nan;status[:,:,j]=1
                    else:
                        l=losses(predict(family,a,yt),ct,family['P']).reshape(nt,reps,3)
                        bad=~np.isfinite(l).all(axis=-1)
                        loss[:,:,j]=l;status[:,:,j]=bad.astype(np.uint8)
                count=int(np.count_nonzero(status));reconstructions+=nt*reps*len(all_candidates)
                if count: failures.append(f'n{N}_g{game}_{scenario}: {count} failures')
                np.savez_compressed(out/f'n{N}_g{game}_{scenario}_loss.npz',loss=loss,status=status)
                # Equality of projected errors at shared ridge values.
                err=max(np.max(abs(loss[:,:,8+1,1]-loss[:,:,8+2,1])),np.max(abs(loss[:,:,8+3,1]-loss[:,:,8+4,1])))
                if not np.isfinite(err) or err>1e-9: failures.append(scenario+': loss projection identity')
                diagnostics.append(dict(N=N,game=game,scenario=scenario,failures=count,identity_error=float(err)))
                if sum(f.stat().st_size for f in BASE.rglob('*') if f.is_file())>180_000_000:
                    raise RuntimeError('Reserved 20 MB for documents; output stop')
    elapsed=time.process_time()-start
    used=sum(f.stat().st_size for f in out.rglob('*') if f.is_file())
    factor=(15*10*400*4)/(1*4*24*2)
    report=dict(success=not failures,failures=failures,cpu_seconds=elapsed,wall_seconds=time.monotonic()-wall,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        bytes=used,reconstructions=reconstructions,max_operator_projection_error=max_identity,
        diagnostics=diagnostics)
    if pilot:
        report.update(estimated_final_cpu_with_10x_margin=elapsed*factor*10,
                      estimated_final_bytes_raw=15*10*400*4*(13*3*8+13)+25_000_000)
        report['feasibility_pass']=bool(not failures and report['estimated_final_cpu_with_10x_margin']<14000 and
            report['estimated_final_bytes_raw']<180_000_000 and report['peak_rss_bytes']<2*1024**3)
    dump(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in ['diagnostics']},indent=2))
    if failures: raise SystemExit(1)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');args=ap.parse_args();run(args.pilot)

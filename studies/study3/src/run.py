"""Frozen study 3 campaign; no study 1/2 output is written."""
import argparse,hashlib,json,time,resource,sys
from pathlib import Path
import numpy as np
import scipy
from laplace_learning.direct import fd_matrix
from studies.study3.src.data import *
from studies.study3.src.model import fit_candidates,select,predict,METHODS,FIXED
BASE=Path(__file__).resolve().parents[1]
ROOTDIR=BASE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def meta(c):
 return {k:(v if not isinstance(v,float) or np.isfinite(v) else None) for k,v in c.items() if k not in ['R','centre']}
def run(pilot=False):
 start=time.process_time();wall=time.monotonic()
 frozen=json.loads((BASE/'FREEZE.json').read_text())
 for p,h in frozen['sha256'].items():assert sha(BASE/p)==h
 tests=json.loads((BASE/'results/tests.json').read_text());assert tests['success']
 for p,h in tests['source_hashes'].items():assert sha(BASE/p)==h
 out=BASE/'results'/('pilot' if pilot else 'final')
 if out.exists():raise RuntimeError('Existing campaign preserved; refusing overwrite')
 if not pilot:
  prior=json.loads((BASE/'results/pilot/report.json').read_text());assert prior['success'] and prior['feasibility_pass']
 out.mkdir(parents=True)
 ns=[32] if pilot else [32,64,256];games=list(range(1 if pilot else 5))
 nv,nt,reps=(24,24,2) if pilot else (256,400,4)
 scenarios=[SCENARIOS[i] for i in [0,1,6,9]] if pilot else SCENARIOS
 offset=100 if pilot else 0
 cfg=dict(pilot=pilot,root_seed=ROOT,sizes=ns,games=games,validation=nv,test=nt,repetitions=reps,
  scenarios=scenarios,methods=METHODS,fixed_methods=FIXED,loss_axes=['L2_squared','modes_ge8_squared','nominal_residual_squared'],
  protocol_sha256=sha(BASE/'PROTOCOL.md'),freeze_commit='8c23fccb9e00d053af71fa9349c1b761c30ca4b4',
  source_sha256={str(p.relative_to(BASE)):sha(p) for p in (BASE/'src').glob('*.py')},
  reused_source_sha256={p:sha(ROOTDIR/p) for p in ['laplace_learning/direct.py','studies/study2/src/model.py']},
  numpy=np.__version__,scipy=scipy.__version__,python=sys.version,
  candidate_counts={'S':140,'M0':20,'Diag':20,'E':20,'H1':20,'H1m':20},status_codes={'success':0,'invalid_estimator':1,'nonfinite_prediction_or_loss':2})
 dump(out/'configuration.json',cfg)
 A_full=fd_matrix(64);failures=[];records=[];fitcount=0;recons=0;invalid_candidates=0
 truth={2.5:profiles(30+offset,0,nt),1.5:profiles(31+offset,0,nt,p=1.5)}
 eps={p:noises(part+offset,0,nt,reps) for p,part in [(2.5,30),(1.5,31)]}
 np.savez_compressed(out/'truth.npz',nominal=truth[2.5],rough=truth[1.5])
 for game in games:
  tr=profiles(10+offset,game,max(ns));tn=noises(10+offset,game,max(ns),1)
  val=profiles(20+offset,game,nv);vn=noises(20+offset,game,nv,2)
  for N in ns:
   fitted={}
   for scenario,m,sigma,p,H in scenarios:
    idx=mask(m);key=(m,sigma);A=A_full[idx]/8
    if key not in fitted:
     yt=observations(tr[:N],tn[:N],idx,sigma)[:,0]
     family=fit_candidates(A,tr[:N,:64],yt,sigma**2/64)
     yv=observations(val,vn,idx,sigma).reshape(-1,len(idx))
     selected,fixed,scores=select(family,yv,np.repeat(val,2,axis=0))
     invalid_candidates+=sum(not x['valid'] for x in scores);fitcount+=1
     fitted[key]=(family,selected,fixed)
     label=f'n{N}_g{game}_{m}_{sigma}'
     dump(out/(label+'_selection.json'),dict(selected={k:meta(v) for k,v in selected.items()},diagnostics=[meta(a) for a in fixed],candidates=[meta(a) for a in scores],
       trace_Cemp=float(np.trace(family['Cemp'])),s=float(family['s']),centre_difference_norm=float(np.linalg.norm(family['my']-A@family['mc']))))
     np.savez_compressed(out/(label+'_models.npz'),R=np.stack([selected[k]['R'] for k in METHODS]+[a['R'] for a in fixed]),
       centres=np.stack([selected[k]['centre'] for k in METHODS]+[a['centre'] for a in fixed]),mc=family['mc'],my=family['my'],A=A,Cemp=family['Cemp'],CH1=family['CH1'],Q=family['Q'],s=family['s'])
    family,selected,fixed=fitted[key]
    y=observations(truth[p],eps[p],idx,sigma,H).reshape(-1,len(idx));c=np.repeat(truth[p],reps,axis=0)
    candidates=[selected[k] for k in METHODS]+fixed
    loss=np.empty((nt,reps,len(candidates),3));status=np.zeros((nt,reps,len(candidates)),dtype=np.uint8)
    for j,a in enumerate(candidates):
     if not a['valid']:loss[:,:,j]=np.nan;status[:,:,j]=1;continue
     pred=predict(family,a,y);v=losses(pred,c,y,A).reshape(nt,reps,3)
     loss[:,:,j]=v;status[:,:,j]=np.where(np.isfinite(v).all(axis=-1),0,2)
    bad=int(np.count_nonzero(status));recons+=nt*reps*len(candidates)
    if bad:failures.append(f'n{N}_g{game}_{scenario}: {bad} entries')
    np.savez_compressed(out/f'n{N}_g{game}_{scenario}_loss.npz',loss=loss,status=status)
    records.append(dict(N=N,game=game,scenario=scenario,failures=bad))
    if sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())>160_000_000:raise RuntimeError('Output reserve guard')
 elapsed=time.process_time()-start;size=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
 report=dict(success=not failures,failures=failures,invalid_candidates=invalid_candidates,conditions=len(records),fits=fitcount,reconstructions=recons,
  cpu_seconds=elapsed,wall_seconds=time.monotonic()-wall,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,bytes=size,records=records)
 if pilot:
  factor=(15*10*400*4)/(4*24*2)
  report.update(estimated_final_cpu_with_10x_margin=elapsed*factor*10,estimated_final_output_bytes=15*10*400*4*(8*3*8+8)+35_000_000)
  report['feasibility_pass']=bool(not failures and not invalid_candidates and report['estimated_final_cpu_with_10x_margin']<14000 and report['estimated_final_output_bytes']<180_000_000 and report['peak_rss_bytes']<2*1024**3)
 dump(out/'report.json',report);print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
 if failures:raise SystemExit(1)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');run(ap.parse_args().pilot)

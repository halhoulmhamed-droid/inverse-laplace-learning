"""Measured global study-3 CPU ledger outside Git; hard process resource limits."""
import json,os,resource,subprocess,sys,time
from pathlib import Path
private=Path('/workspace/inverse-laplace-private/study3');private.mkdir(parents=True,exist_ok=True)
ledger=private/'resources.json'
entries=json.loads(ledger.read_text()) if ledger.exists() else []
used=sum(e['cpu_seconds'] for e in entries)
remaining=int(14400-used-30)
if remaining<=0: raise SystemExit('CPU budget exhausted')
if any(e.get('running') for e in entries): raise SystemExit('Unclosed prior reservation: inspect privately')
entry=dict(command=sys.argv[1:],running=True,reserved_cpu=remaining)
entries.append(entry);ledger.write_text(json.dumps(entries,indent=2))
def limits():
 resource.setrlimit(resource.RLIMIT_CPU,(remaining,remaining))
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR='/tmp/study3-mpl')
a=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
log=private/f'command-{len(entries):02d}.log'
with log.open('w') as f:
 p=subprocess.run(sys.argv[1:],env=env,preexec_fn=limits,stdout=f,stderr=subprocess.STDOUT)
b=resource.getrusage(resource.RUSAGE_CHILDREN)
entry.update(running=False,exit_code=p.returncode,cpu_seconds=b.ru_utime+b.ru_stime-a.ru_utime-a.ru_stime,
             wall_seconds=time.monotonic()-start,peak_rss_bytes=b.ru_maxrss*1024)
ledger.write_text(json.dumps(entries,indent=2)+'\n')
print(log.read_text()[-5000:]);print(json.dumps(entry))
raise SystemExit(p.returncode)

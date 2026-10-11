"""Independent scalar recomputation from seeds and stored affine maps; not a new campaign."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parents[1];root=BASE.parents[1]
result={};cfg=json.loads((BASE/'results/final/configuration.json').read_text())
for p,h in cfg['sources'].items():assert hashlib.sha256((BASE/p).read_bytes()).hexdigest()==h
assert hashlib.sha256((root/'laplace_learning/direct.py').read_bytes()).hexdigest()==cfg['reused_direct_sha256']
fr=json.loads((BASE/'FREEZE.json').read_text())
for p,h in fr['sha256'].items():assert hashlib.sha256((BASE/p).read_bytes()).hexdigest()==h
changes=subprocess.check_output(['git','diff',fr['baseline_commit'],'--name-only','--','. ',':(exclude)studies/study2/**'],cwd=root,text=True)
# Exact baseline tree comparison, including old PDFs and results.
files=subprocess.check_output(['git','ls-tree','-r','--name-only',fr['baseline_commit']],cwd=root,text=True).splitlines()
for p in files:
 expected=subprocess.check_output(['git','show',fr['baseline_commit']+':'+p],cwd=root)
 assert (root/p).read_bytes()==expected,p
result['study1_files_byte_identical']=len(files)
checks=0;max_error=0
for scenario,p,height,partition in [('primary',2.5,.25,30),('rough',1.5,.25,31),('height10',2.5,.275,30)]:
 z=np.load(BASE/'results/final/n256_g0_half_0.01_models.npz')
 maps=np.concatenate([z['R'],z['fixed_R']]);V=z['V'];P=V@V.T
 raw=np.load(BASE/f'results/final/n256_g0_{scenario}_loss.npz')['loss']
 for obj in [0,399]:
  weights=np.array([(1+k*k)**(-p/2) for k in range(128)]);weights/=np.linalg.norm(weights)
  c=np.random.Generator(np.random.PCG64(np.random.SeedSequence([2026101102,partition,0,obj,0,0]))).normal(size=128)*weights
  for rep in [0,3]:
   eps=np.random.Generator(np.random.PCG64(np.random.SeedSequence([2026101102,partition,0,obj,1,rep]))).normal(size=64)
   # Integrate cosine means directly with endpoint sines, not the production sinc formula.
   B=np.ones((32,128))
   for j in range(32):
    for k in range(1,128):
     B[j,k]=np.sqrt(2)*64*(np.sin(k*np.pi*(j+1)/64)-np.sin(k*np.pi*j/64))/(k*np.pi*np.cosh(k*np.pi*height))
   y=(B@c+.01*eps[:32])/8
   for j,R in enumerate(maps):
    pred=z['mc']+R@(y-z['my']);e=pred-c[:64]
    full=np.sum((np.concatenate([pred,np.zeros(64)])-c)**2)
    parts=np.array([full,np.dot(P@e,P@e),np.dot(e-P@e,e-P@e)+np.dot(c[64:],c[64:])])
    error=float(np.max(abs(parts-raw[obj,rep,j])));max_error=max(max_error,error);checks+=3
    np.testing.assert_allclose(parts,raw[obj,rep,j],rtol=2e-10,atol=1e-11)
result.update(independent_scalar_checks=checks,max_absolute_difference=max_error)
report=json.loads((BASE/'results/analysis/report.json').read_text())
for p,h in report['raw_sha256'].items():assert hashlib.sha256((BASE/'results/final'/p).read_bytes()).hexdigest()==h
result['raw_loss_files_sha256_verified']=len(report['raw_sha256'])
result['success']=True
(BASE/'results/verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

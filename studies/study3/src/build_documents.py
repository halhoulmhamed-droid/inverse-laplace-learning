"""Compile study 3 only; operational logs outside the public repository."""
from pathlib import Path
import subprocess,shutil,json
BASE=Path(__file__).resolve().parents[1]
OUT=Path('/workspace/inverse-laplace-private/study3/latex');OUT.mkdir(parents=True,exist_ok=True)
for name in ['report','slides']:
 for iteration in [1,2]:
  p=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory='+str(OUT),name+'.tex'],cwd=BASE/'docs',capture_output=True,text=True,errors='replace')
  (OUT/f'{name}-pass{iteration}.txt').write_text(p.stdout+p.stderr)
  if p.returncode:print(p.stdout[-4000:]);raise SystemExit(p.returncode)
 shutil.copy2(OUT/(name+'.pdf'),BASE/'docs'/(name+'.pdf'))
print(json.dumps({'success':True,'documents':['docs/report.pdf','docs/slides.pdf']}))

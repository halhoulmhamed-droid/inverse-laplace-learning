"""Compile the French report and slides, preserving logs in an ignored build dir."""
from pathlib import Path
import json
import resource
import shutil
import subprocess

root=Path(__file__).resolve().parents[1]
out=root/'build'/'documents';out.mkdir(parents=True,exist_ok=True)
start=resource.getrusage(resource.RUSAGE_CHILDREN)
for stem in ['report','slides']:
    for pass_number in range(2):
        result=subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',
            '-output-directory='+str(out),stem+'.tex'],cwd=root/'docs',capture_output=True,text=True,errors="replace")
        (out/f'{stem}.pass{pass_number+1}.txt').write_text(result.stdout+result.stderr)
        if result.returncode:
            print(result.stdout[-4000:]);raise SystemExit(result.returncode)
    shutil.copy2(out/(stem+'.pdf'),root/'docs'/(stem+'.pdf'))
end=resource.getrusage(resource.RUSAGE_CHILDREN)
report={'success':True,'documents':['docs/report.pdf','docs/slides.pdf'],
        'cpu_seconds':end.ru_utime+end.ru_stime-start.ru_utime-start.ru_stime,
        'compiler':subprocess.check_output(['pdflatex','--version'],text=True).splitlines()[0]}
(out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

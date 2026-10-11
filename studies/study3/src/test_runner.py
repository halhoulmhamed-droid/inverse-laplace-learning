"""Run targeted tests and bind the prerequisite report to current source bytes."""
import unittest,json,hashlib,time
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
t=time.process_time()
suite=unittest.defaultTestLoader.discover(str(BASE/'tests'),pattern='test_*.py')
r=unittest.TextTestRunner(verbosity=2).run(suite)
report=dict(success=r.wasSuccessful(),tests_run=r.testsRun,failures=len(r.failures),errors=len(r.errors),cpu_seconds=time.process_time()-t,
 source_hashes={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for folder in ['src','tests'] for p in (BASE/folder).glob('*.py')})
(BASE/'results/tests.json').write_text(json.dumps(report,indent=2)+'\n')
raise SystemExit(not r.wasSuccessful())

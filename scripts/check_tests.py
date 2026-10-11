"""Execute independent unit checks and record the exact scientific source hashes."""
import hashlib
import json
from pathlib import Path
import resource
import unittest
import argparse

root=Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0,str(root))
parser=argparse.ArgumentParser()
parser.add_argument('--output',default=str(root/'results'/'tests.json'))
args=parser.parse_args()
start=resource.getrusage(resource.RUSAGE_SELF)
suite=unittest.defaultTestLoader.discover(str(root/'tests'))
result=unittest.TextTestRunner(verbosity=2).run(suite)
end=resource.getrusage(resource.RUSAGE_SELF)
evidence=dict(success=result.wasSuccessful(),tests_run=result.testsRun,
 cpu_seconds=end.ru_utime+end.ru_stime-start.ru_utime-start.ru_stime,
 source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'laplace_learning').glob('*.py'))},
 test_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((root/'tests').glob('test_*.py'))})
Path(args.output).write_text(json.dumps(evidence,indent=2)+'\n')
raise SystemExit(not result.wasSuccessful())

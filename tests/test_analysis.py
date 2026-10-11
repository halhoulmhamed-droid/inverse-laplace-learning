import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from scripts.analyze import clean, run, select_diagnostics, summarize


class AnalysisTests(unittest.TestCase):
    def test_extended_summary(self):
        self.assertEqual(summarize([1, 2, np.inf])['q95'], np.inf)
        self.assertEqual(clean({'a': np.inf, 'b': None}), {'a': '+inf', 'b': None})

    def test_diagnostic_ties_and_failure(self):
        e = np.full((20, 2), 2.)
        r = np.full((20, 2), 3.)
        flags = np.zeros((20, 2, 11), dtype='uint8')
        flags[3, 1, 4] = 1
        d = select_diagnostics(e, r, flags)
        self.assertEqual([d[k] for k in ('median', 'q95', 'worst')], [0, 0, 0])
        self.assertEqual(d['low_residual_high_error'], 0)
        self.assertEqual(d['failure'], [3, 1])
        e[0] = np.inf
        r[0] = np.inf
        self.assertNotEqual(select_diagnostics(e, r, flags)['low_residual_high_error'], 0)

    def test_no_analysis_of_incomplete_campaign(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)
            (p/'report.json').write_text(json.dumps({'verdict': 'PASS', 'complete': False}))
            with self.assertRaisesRegex(ValueError, 'complete'):
                run(p, p/'out')
            self.assertFalse((p/'out').exists())

if __name__ == '__main__':
    unittest.main()

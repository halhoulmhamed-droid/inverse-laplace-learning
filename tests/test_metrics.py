import unittest

import numpy as np

from laplace_learning.metrics import (
    bootstrap_ratios, metrics_for_predictions, ordered_quantile,
    paired_ratio_summary, relative_errors,
)


class MetricsTests(unittest.TestCase):
    def test_ordered_quantile_keeps_infinity(self):
        self.assertEqual(ordered_quantile([1., 2., 3., np.inf], .5), 2.)
        self.assertEqual(ordered_quantile([1., np.inf], .95), np.inf)
        with self.assertRaises(ValueError):
            ordered_quantile([1., np.nan], .5)

    def test_ratio_has_no_false_win_against_failed_baseline(self):
        for baseline in (0., np.inf, np.nan):
            self.assertFalse(paired_ratio_summary(1., baseline)["defined"])
        self.assertEqual(paired_ratio_summary(np.inf, 2.)["ratio"], np.inf)
        np.testing.assert_equal(relative_errors([[1., np.inf], [2., 4.]],
                                               [0., 2.], .1),
                                [[10., np.inf], [1., 2.]])

    def test_visible_kernel_and_unrepresented_tail(self):
        out = metrics_for_predictions([[1., 3.]], [[0., 0., 4.]],
                                      [[1., 0.]], [[2.]], [[5.]], [[1.], [0.]])
        self.assertAlmostEqual(out["error_abs"][0], np.sqrt(26.))
        self.assertEqual(out["visible"][0], 1.)
        self.assertEqual(out["null"][0], 3.)
        self.assertAlmostEqual(out["error_abs"][0] ** 2,
                               out["visible"][0] ** 2 + out["null"][0] ** 2 + 4. ** 2)
        self.assertEqual(out["residual_nominal"][0], 1.)
        self.assertEqual(out["residual_independent"][0], 4.)
        self.assertEqual(out["status"][0], 0)

    def test_scaled_norm_and_distinct_failures(self):
        out = metrics_for_predictions([[1.e200], [np.nan], [1.e308]],
                                      [[0.], [0.], [-1.e308]], [[0.]],
                                      [[0.], [0.], [0.]], [[0.], [0.], [0.]], [[1.]])
        self.assertEqual(out["error_abs"][0], 1.e200)
        np.testing.assert_equal(out["status"], [0, 1, 2])
        self.assertTrue(np.isinf(out["error_abs"][1:]).all())
        with self.assertRaisesRegex(ValueError, "incomplete"):
            metrics_for_predictions(None, [[0.]], [[1.]], [[0.]], [[0.]], [[1.]])
        with self.assertRaisesRegex(ValueError, "corrupt"):
            metrics_for_predictions([[0.]], [[np.nan]], [[1.]], [[0.]], [[0.]], [[1.]])

    def test_cluster_draw_matches_explicit_seed_and_is_reproducible(self):
        base = np.array([[1., 2.], [4., 8.], [16., 32.]])
        capped = base * np.array([[.2], [.5], [.8]])
        result = bootstrap_ratios(capped, base, capped, base, 1, n_boot=20)
        again = bootstrap_ratios(capped, base, capped, base, 1, n_boot=20)
        np.testing.assert_equal(result["draws"], again["draws"])
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(
            [20261011, 5, 1, 0, 2, 0])))
        ids = rng.integers(0, 3, size=3)
        expected = ordered_quantile(capped[ids], .95) / ordered_quantile(base[ids], .95)
        self.assertEqual(result["draws"][0, 0], expected)

    def test_invalid_bootstrap_draws_are_retained(self):
        base = np.array([[0., 0.], [1., 1.]])
        result = bootstrap_ratios(base * .5, base, base * .5, base, 0, n_boot=100)
        self.assertGreater(result["invalid_q95"], 0)
        self.assertIsNone(result["upper_q95"])
        self.assertEqual(result["verdict"], "NON_CONCLUSIVE")
        self.assertEqual(result["draws"].shape, (100, 2))
        failed = bootstrap_ratios(np.full((2, 2), np.inf), np.ones((2, 2)),
                                  np.ones((2, 2)), np.ones((2, 2)), 0, n_boot=10)
        self.assertEqual(failed["upper_q95"], np.inf)
        self.assertEqual(failed["verdict"], "NOT_SATISFIED")

    def test_infinite_point_estimate_cannot_pass_finite_small_bootstrap(self):
        # With this frozen seed the single draw selects object 1 twice.
        # It misses the failed object, but the point q95 is still infinite.
        result = bootstrap_ratios([[np.inf], [.5]], [[1.], [1.]],
                                  [[.5], [.5]], [[1.], [1.]], 0, n_boot=1)
        self.assertEqual(result["q95_ratio"], np.inf)
        self.assertEqual(result["upper_q95"], .5)
        self.assertEqual(result["verdict"], "NOT_SATISFIED")


if __name__ == "__main__":
    unittest.main()

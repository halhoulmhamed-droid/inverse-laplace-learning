"""Independent algebraic checks; unittest needs no additional dependencies."""
import unittest
import numpy as np
from numpy.testing import assert_allclose
from laplace_learning.filters import fit_candidates, select_candidates, predict, METHODS


class FilterTests(unittest.TestCase):
    def setUp(self):
        self.c = np.array([[-1., -2.], [1., 2.], [-2., 1.], [2., -1.]])
        self.A = np.array([[1., .4], [0., .3]])
        self.d = self.c @ self.A.T + np.array([3., -2.])
        self.candidates = fit_candidates(self.A, self.c, self.d, .1)

    def one(self, method, index=0):
        return next(x for x in self.candidates if x['method'] == method and x['grid_index'] == index)

    def test_counts_and_grid(self):
        self.assertEqual(set(x['method'] for x in self.candidates), set(METHODS))
        self.assertEqual(len(self.candidates), 126)  # 20*6 + 2 ranks + 4 fixed
        candidates = fit_candidates(np.eye(20), np.eye(20), np.eye(20), .01)
        self.assertEqual(len(candidates), 144)

    def test_noncommuting_h1_normal_equations(self):
        est = self.one('tikhonov_h1', 12)
        P = np.diag(np.array([1., 1. + np.pi**2]) / (1. + np.pi**2))
        self.assertGreater(np.linalg.norm(P @ self.A.T @ self.A - self.A.T @ self.A @ P), .01)
        assert_allclose((self.A.T @ self.A + est['parameter'] * P) @ est['R'], self.A.T, atol=1e-13)
        self.assertAlmostEqual(est['operator_norm'], np.linalg.svd(est['R'])[1][0])

    def test_affine_offset(self):
        est = self.one('learned_signed')
        assert_allclose(predict(est, est['mu_d']), est['mu_c'])
        truth = np.array([1.2, -.7])
        data = self.A @ truth + np.array([.3, -.2])
        offset = self.A @ est['mu_c'] - est['mu_d']
        expected = (est['R'] @ self.A - np.eye(2)) @ (truth - est['mu_c']) + est['R'] @ (data - self.A @ truth + offset)
        assert_allclose(predict(est, data) - truth, expected)

    def test_negative_covariance_and_weak_denominator(self):
        c = np.array([[-1.], [1.]])
        candidates = fit_candidates(np.ones((1, 1)), c, -c + 5, 0.)
        by = {x['method']: x for x in candidates}
        assert_allclose(by['learned_signed']['R'], [[-1.]])
        assert_allclose(by['learned_positive']['R'], [[0.]])
        nearzero = fit_candidates(np.ones((1, 1)), c, 1e-20 * c, 0.)
        fixed = next(x for x in nearzero if x['method'] == 'learned_signed')
        assert_allclose(fixed['R'], [[0.]])
        self.assertEqual(fixed['rejected_variances'], 1)

    def test_wiener_equals_learned_for_uncorrelated_sample(self):
        # Exact zero means and sample cross covariance; p=v=4/3, sigma=.5.
        c = np.array([[-1.], [-1.], [1.], [1.]])
        noise = np.array([[-1.], [1.], [-1.], [1.]])
        candidates = fit_candidates(np.array([[.5]]), c, .5*c + noise, 4/3)
        by = {x['method']: x for x in candidates}
        assert_allclose(by['wiener']['R'], [[.4]])
        assert_allclose(by['learned_signed']['R'], by['wiener']['R'])
        # No signal/no noise: zero gain is a defined convention.
        null = fit_candidates(np.ones((1, 1)), np.ones((3, 1)), np.ones((3, 1)), 0.)
        self.assertTrue(all(x['valid'] for x in null))

    def test_lmmse_predicts_correlated_null_component(self):
        c = np.array([[-1., -2.], [1., 2.]]) + [3., 7.]
        d = c[:, :1] + 4
        candidates = fit_candidates(np.array([[1., 0.]]), c, d, 0.)
        est = next(x for x in candidates if x['method'] == 'lmmse' and x['grid_index'] == 0)
        expected = np.array([[1.], [2.]]) / (1 + est['parameter'])
        assert_allclose(est['R'], expected)
        assert_allclose(predict(est, np.array([[7.]])), [[3., 7.]])

    def test_full_tikh_keeps_numerically_discarded_direction(self):
        A = np.diag([1., 1e-13])
        candidate = fit_candidates(A, self.c, self.c @ A.T, .1)[0]
        self.assertEqual(candidate['numerical_rank'], 1)
        self.assertGreater(candidate['R'][1, 1], 0.)
        assert_allclose(candidate['R'][1, 1], 1e-13/(1e-26 + candidate['parameter']))

    def test_validation_nonfinite_and_ties(self):
        def item(R, index):
            return dict(method='example', parameter=index, grid_index=index,
                        R=np.array([[R]]), mu_c=np.zeros(1), mu_d=np.zeros(1),
                        operator_norm=abs(R), valid=True)
        candidates = [item(2., 0), item(1., 2), item(1., 1), item(np.nan, 3)]
        selected = select_candidates(candidates, np.ones((1, 1)), np.zeros((2, 1)), np.zeros((2, 1)))
        self.assertEqual(selected['example']['grid_index'], 1)
        invalid = select_candidates(candidates, np.ones((1, 1)), np.full((2, 1), np.nan), np.zeros((2, 1)))['example']
        self.assertFalse(invalid['valid'])
        self.assertEqual(invalid['failure'], 'all_candidates_invalid')
        self.assertTrue(np.isnan(predict(invalid, np.zeros((1, 1)))).all())

    def test_overflow_is_not_repaired_by_cap(self):
        c = np.array([[-1e200], [1e200]])
        candidates = fit_candidates(np.ones((1, 1)), c, c, .1)
        learned = [x for x in candidates if x['method'].startswith('learned')]
        self.assertTrue(all(not x['valid'] for x in learned))
        self.assertTrue(all(not x['valid'] for x in candidates if x['method'].startswith('wiener')))

    def test_validation_diagnostics(self):
        selected, diagnostics = select_candidates(
            self.candidates, self.A, self.d, self.d, return_diagnostics=True)
        self.assertEqual(len(diagnostics), len(self.candidates))
        self.assertEqual(set(selected), set(METHODS))
        self.assertTrue(all('R' not in item for item in diagnostics))
        for estimator in selected.values():
            metadata = next(item for item in diagnostics
                            if item['method'] == estimator['method']
                            and item['grid_index'] == estimator['grid_index'])
            self.assertEqual(metadata['score'], estimator['score'])
            self.assertEqual(metadata['valid'], estimator['valid'])

    def test_zero_rank_and_input_validation(self):
        candidates = fit_candidates(np.zeros((2, 2)), self.c, np.ones((4, 2)), .1)
        tsvd = [x for x in candidates if x['method'] == 'tsvd']
        self.assertEqual(len(tsvd), 1)
        assert_allclose(tsvd[0]['R'], np.zeros((2, 2)))
        with self.assertRaises(ValueError):
            fit_candidates(self.A, self.c[:1], self.d[:1], .1)


if __name__ == '__main__':
    unittest.main()

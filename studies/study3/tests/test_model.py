"""Independent identities and selection checks for the frozen study 3 families."""
import unittest
from unittest.mock import patch

import numpy as np

from studies.study3.src import model
from studies.study2.src.model import covariance_family as study2_family
from studies.study2.src.model import predict as study2_predict


class CovarianceModelAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(73003)
        cls.A = rng.normal(size=(5, 7))
        cls.c = rng.normal(size=(36, 7)) + np.linspace(-.3, .8, 7)
        cls.y = cls.c @ cls.A.T + .2 * rng.normal(size=(36, 5)) + .3
        cls.nu = .04
        cls.family = model.fit_candidates(cls.A, cls.c, cls.y, cls.nu)

    def candidate(self, method, index=8, family=None):
        source = self.family if family is None else family
        return next(c for c in source['candidates']
                    if c['method'] == method and c['index'] == index)

    def test_covariance_trace_psd_and_target_scale(self):
        f = self.family
        np.testing.assert_allclose(f['Cemp'], np.cov(self.c, rowvar=False), atol=1e-14)
        energy = np.trace(f['Cemp'])
        q_inverse = np.linalg.inv(f['Q'])
        np.testing.assert_allclose(f['CH1'], energy * q_inverse / np.trace(q_inverse))
        scaled_inverse = np.linalg.inv(13 * f['Q'])
        np.testing.assert_allclose(f['CH1'], energy * scaled_inverse / np.trace(scaled_inverse))
        for gamma in model.GAMMAS:
            C = (1 - gamma) * f['Cemp'] + gamma * f['CH1']
            np.testing.assert_allclose(C, C.T, atol=1e-14)
            self.assertAlmostEqual(np.trace(C), energy, places=12)
            self.assertGreater(np.linalg.eigvalsh(C).min(), -1e-12)

    def test_gamma_one_matches_exact_h1_parameter_and_centre(self):
        f = self.family
        for j in [0, 8, 19]:
            beta = model.ALPHAS[j] * self.nu
            eta = beta / f['s']
            exact = np.linalg.solve(self.A.T @ self.A + eta * f['Q'], self.A.T)
            np.testing.assert_allclose(self.candidate('H1m', j)['R'], exact,
                                       atol=2e-10, rtol=2e-10)
            endpoint = self.candidate('S', 6 * len(model.ALPHAS) + j)
            np.testing.assert_array_equal(endpoint['R'], self.candidate('H1m', j)['R'])
            np.testing.assert_allclose(endpoint['centre'], self.A @ f['mc'])
            qn = f['Q'] / f['Q'].max()
            historical_a = beta * f['Q'].max() / f['s']
            exact_normalised = np.linalg.solve(self.A.T @ self.A + historical_a * qn,
                                               self.A.T)
            np.testing.assert_allclose(endpoint['R'], exact_normalised,
                                       atol=2e-10, rtol=2e-10)

    def test_gamma_zero_rank_deficiency_and_observation_residual(self):
        f = model.fit_candidates(self.A, self.c[:3], self.y[:3], self.nu)
        self.assertLessEqual(np.linalg.matrix_rank(f['Cemp']), 2)
        _, _, vh = np.linalg.svd(self.c[:3] - self.c[:3].mean(0), full_matrices=True)
        null = vh[2:].T
        for j in [0, 8, 19]:
            candidate = self.candidate('M0', j, f)
            self.assertTrue(candidate['valid'])
            R = candidate['R']
            system = self.A @ f['Cemp'] @ self.A.T + model.ALPHAS[j] * self.nu * np.eye(5)
            self.assertGreater(np.linalg.eigvalsh(system).min(), 0)
            np.testing.assert_allclose(R @ system, f['Cemp'] @ self.A.T,
                                       atol=1e-10, rtol=1e-10)
            np.testing.assert_allclose(null.T @ R, 0, atol=1e-10)
            np.testing.assert_array_equal(R, self.candidate('S', j, f)['R'])

    def test_full_rank_plugin_solves_variational_normal_equation(self):
        f = self.family
        gamma_index, alpha_index = 3, 9
        gamma = model.GAMMAS[gamma_index]
        candidate = self.candidate('S', gamma_index * len(model.ALPHAS) + alpha_index)
        C = (1 - gamma) * f['Cemp'] + gamma * f['CH1']
        precision = np.linalg.inv(C)
        beta = model.ALPHAS[alpha_index] * self.nu
        y = self.y[:4]
        increment = model.predict(f, candidate, y) - f['mc']
        derivative = (increment @ self.A.T - (y - self.A @ f['mc'])) @ self.A
        derivative += beta * increment @ precision
        np.testing.assert_allclose(derivative, 0, atol=1e-11)

    def test_diagonal_empirical_covariance_is_same_family(self):
        rng = np.random.default_rng(517)
        orthogonal = np.linalg.qr(np.column_stack([np.ones(20), rng.normal(size=(20, 7))]))[0][:, 1:]
        c = orthogonal * np.linspace(.2, 2., 7)
        y = c @ self.A.T
        f = model.fit_candidates(self.A, c, y, self.nu)
        np.testing.assert_allclose(f['Cemp'], np.diag(np.diag(f['Cemp'])), atol=1e-16)
        for j in [0, 8, 19]:
            np.testing.assert_allclose(self.candidate('M0', j, f)['R'],
                                       self.candidate('Diag', j, f)['R'],
                                       atol=1e-11, rtol=1e-11)

    def test_empirical_and_historical_h1_match_study2_exactly(self):
        old = study2_family(self.A, self.c, self.y, rank=3)
        for name, old_name in [('E', 'L'), ('H1', 'H1')]:
            for j in range(20):
                current = self.candidate(name, j)
                previous = next(c for c in old['candidates']
                                if c['method'] == old_name and c['index'] == j)
                self.assertEqual(current['parameter'], previous['parameter'])
                np.testing.assert_allclose(current['R'], previous['R'], atol=1e-13, rtol=1e-13)
                np.testing.assert_array_equal(current['centre'], old['my'])
                np.testing.assert_allclose(model.predict(self.family, current, self.y[:3]),
                                           study2_predict(old, previous, self.y[:3]),
                                           atol=1e-13, rtol=1e-13)

    def test_empirical_and_plugin_differ_for_finite_noisy_pairs(self):
        empirical = self.candidate('E')
        plugin = self.candidate('M0')
        self.assertGreater(np.linalg.norm(empirical['centre'] - plugin['centre']), .1)
        self.assertGreater(np.linalg.norm(empirical['R'] - plugin['R']), 1e-4)
        reconstruction = model.predict(self.family, plugin, self.y[:3])
        erroneous_centre = self.family['mc'] + (self.y[:3] - self.family['my']) @ plugin['R'].T
        self.assertGreater(np.linalg.norm(reconstruction - erroneous_centre), 1e-4)

    def test_noiseless_empirical_plugin_identity_with_converted_ridge(self):
        y = self.c @ self.A.T
        ty = np.trace(np.cov(y, rowvar=False)) / y.shape[1]
        empirical_index, alpha_index = 9, 8
        noise_variance = (model.EMPIRICAL_LAMBDAS[empirical_index] * ty
                          / model.ALPHAS[alpha_index])
        f = model.fit_candidates(self.A, self.c, y, noise_variance)
        empirical = self.candidate('E', empirical_index, f)
        plugin = self.candidate('M0', alpha_index, f)
        np.testing.assert_allclose(empirical['centre'], plugin['centre'], atol=1e-14)
        np.testing.assert_allclose(empirical['R'], plugin['R'], atol=1e-11, rtol=1e-11)

    def test_validation_selection_budgets_tail_and_matched_parameters(self):
        truth = np.column_stack([self.c[:9], np.ones((9, 4)) * .6])
        selected, fixed, scores = model.select(self.family, self.y[:9], truth)
        self.assertEqual(len(scores), 240)
        self.assertEqual([sum(c['method'] == method for c in scores) for method in model.METHODS],
                         [140, 20, 20, 20, 20, 20])
        self.assertEqual([c['diagnostic'] for c in fixed], model.FIXED)
        self.assertTrue(all(c['alpha_index'] == selected['S']['alpha_index'] for c in fixed))
        for method in model.METHODS:
            chosen = selected[method]
            padded = np.pad(model.predict(self.family, chosen, self.y[:9]), ((0, 0), (0, 4)))
            expected = np.mean(np.sum((padded - truth)**2, axis=1))
            self.assertAlmostEqual(chosen['score'], expected, places=12)
            minimum = min(c['score'] for c in scores if c['method'] == method)
            self.assertLessEqual(chosen['score'] - minimum, 1e-12 * max(1., minimum))
        self.assertNotIn('score', self.family['candidates'][0])

    def test_selection_ties_use_norm_then_index(self):
        candidates = []
        for original in self.family['candidates']:
            current = dict(original, R=np.zeros_like(original['R']), norm=2.)
            if current['index'] in [3, 6]:
                current['norm'] = 1.
            candidates.append(current)
        family = dict(self.family, candidates=candidates)
        selected, _, _ = model.select(family, self.y[:3], self.c[:3])
        self.assertEqual([selected[name]['index'] for name in model.METHODS], [3] * 6)

    def test_solver_failure_is_recorded_without_fallback(self):
        with patch.object(model, '_solve', side_effect=np.linalg.LinAlgError('injected SPD solve failure')):
            family = model.fit_candidates(self.A, self.c, self.y, self.nu)
        self.assertEqual(len(family['candidates']), 240)
        self.assertTrue(all(not c['valid'] and 'injected' in c['failure']
                            for c in family['candidates']))
        selected, fixed, scores = model.select(family, self.y[:3], self.c[:3])
        self.assertTrue(all(not c['valid'] for c in selected.values()))
        self.assertTrue(all(not c['valid'] and np.isnan(c['R']).all() for c in fixed))
        self.assertTrue(all(np.isinf(c['score']) for c in scores))

    def test_invalid_inputs_rejected_explicitly(self):
        for noise in [0., -1., float('inf'), float('nan'), np.array([1., 2.])]:
            with self.assertRaises(ValueError):
                model.fit_candidates(self.A, self.c, self.y, noise)
        for c, y in [(self.c[:1], self.y[:1]), (self.c[:, :3], self.y),
                     (np.full_like(self.c, np.nan), self.y),
                     (np.ones_like(self.c), self.y), (self.c, np.ones_like(self.y))]:
            with self.assertRaises(ValueError):
                model.fit_candidates(self.A, c, y, self.nu)
        with self.assertRaises(ValueError):
            model.select(self.family, self.y[:2], self.c[:3])
        with self.assertRaises(ValueError):
            model.select(self.family, self.y[:2], np.full((2, 7), np.nan))

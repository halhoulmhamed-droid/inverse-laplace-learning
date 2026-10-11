"""Recompute study-3 checks from archives without importing production science.

Run from the repository root through src/bounded.py. No campaign or fitting
selector is rerun. Writes one scientific JSON (results/independent_review.json
by default); operational timing and logs are handled privately by the launcher.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
FINAL = BASE / 'results/final'
SEED = 2026101103


def random(part, game, obj, channel=0, rep=0):
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(
        [SEED, part, game, obj, channel, rep])))


def coefficients(part, game, obj, p=2.5):
    weights = np.array([(1 + k*k)**(-p/2) for k in range(128)])
    weights /= np.linalg.norm(weights)
    return random(part, game, obj).normal(size=128) * weights


def analytic_transfer(height):
    """Interval means from endpoint sine integrals, not production sinc."""
    endpoints = np.arange(65) / 64
    modes = np.arange(1, 128)
    transfer = np.ones((64, 128))
    transfer[:, 1:] = (64 * np.sqrt(2)
        * np.diff(np.sin(np.pi * endpoints[:, None] * modes), axis=0)
        / (np.pi * modes * np.cosh(np.pi * modes * height)))
    return transfer


def plugin(covariance, A, alpha, variance):
    return np.linalg.solve(A @ covariance @ A.T
        + alpha * variance * np.eye(len(A)), A @ covariance).T


def main(output=None):
    output = Path(output) if output is not None else BASE / 'results/independent_review.json'
    if output.exists():
        raise RuntimeError('Existing independent review preserved; choose a private output for reruns')
    cfg = json.loads((FINAL / 'configuration.json').read_text())
    frozen = json.loads((BASE / 'FREEZE.json').read_text())
    analysis = json.loads((BASE / 'results/analysis/report.json').read_text())
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    for p, expected in frozen['sha256'].items():
        assert sha(BASE / p) == expected, p
    for p, expected in cfg['source_sha256'].items():
        assert sha(BASE / p) == expected, p
    for p, expected in cfg['reused_source_sha256'].items():
        assert sha(ROOT / p) == expected, p
    for p, expected in analysis['raw_sha256'].items():
        assert sha(FINAL / p) == expected, p
    methods = cfg['methods'] + cfg['fixed_methods']
    rows = list(csv.DictReader((BASE / 'results/analysis/summary.csv').open()))
    table = {(int(r['N']), int(r['game']), r['scenario'], r['method']): r for r in rows}
    checks, summary_error, failures = 0, 0., 0
    for path in FINAL.glob('*_loss.npz'):
        n, g, scenario = path.stem[:-5].split('_', 2)
        n, g = int(n[1:]), int(g[1:])
        with np.load(path) as archive:
            loss, status = archive['loss'], archive['status']
        assert loss.shape == (400, 4, 8, 3)
        assert np.isfinite(loss).all()
        failures += int(np.count_nonzero(status))
        for j, method in enumerate(methods):
            row = table[n, g, scenario, method]
            observed = [loss[:, :, j, k].mean() for k in range(3)]
            observed.append(np.quantile(loss[:, :, j, 0].mean(axis=1), .95))
            expected = [float(row[k]) for k in ['risk', 'high_mode', 'residual', 'q95']]
            summary_error = max(summary_error, float(np.max(np.abs(np.array(observed)-expected))))
            np.testing.assert_allclose(observed, expected, rtol=1e-14, atol=1e-15)
            checks += 4
    assert failures == analysis['failures'] == 0

    # Shared resampling indices across both methods AND both masks.
    saved_boot = np.load(BASE / 'results/analysis/bootstrap.npz')['ratios']
    decisions, boot_error = [], 0.
    for game in range(5):
        full = np.load(FINAL / f'n64_g{game}_primary_loss.npz')['loss'][:, :, :2, 0].mean(axis=1)
        half = np.load(FINAL / f'n64_g{game}_partial_loss.npz')['loss'][:, :, :2, 0].mean(axis=1)
        rng = random(40, game, 0)
        bootstrap = []
        for _ in range(32):
            indices = rng.integers(0, 400, size=(125, 400))
            x, y = full[indices].sum(axis=1), half[indices].sum(axis=1)
            bootstrap.append(np.column_stack([x[:, 0]/x[:, 1], y[:, 0]/y[:, 1]]))
        bootstrap = np.concatenate(bootstrap)
        upper = np.quantile(bootstrap, .995, axis=0)
        point = np.array([full[:, 0].mean()/full[:, 1].mean(), half[:, 0].mean()/half[:, 1].mean()])
        boot_error = max(boot_error, float(np.max(np.abs(bootstrap-saved_boot[game]))))
        np.testing.assert_allclose(bootstrap, saved_boot[game], rtol=1e-14, atol=1e-15)
        archived = analysis['primary'][game]
        np.testing.assert_allclose(np.r_[point, upper], [archived[k] for k in
            ['primary_ratio', 'protection_ratio', 'primary_upper', 'protection_upper']], rtol=1e-14, atol=1e-15)
        met = bool(upper[0] <= .97 and upper[1] <= 1.03)
        assert met == archived['criteria_met']
        decisions.append(dict(game=game, primary_ratio=float(point[0]),
            protection_ratio=float(point[1]), primary_upper=float(upper[0]),
            protection_upper=float(upper[1]), gain_established=bool(upper[0] <= .97),
            protection_established=bool(upper[1] <= 1.03), joint_criterion_met=met))
    conclusion = 'ESTABLISHED_ON_FROZEN_CONDITIONS' if all(d['joint_criterion_met'] for d in decisions) else 'NOT_ESTABLISHED'
    assert conclusion == analysis['conclusion']

    # Saved matrices: compare covariance formulas, common alpha, H1 endpoint.
    maps_count, zero_count, budget_checks = 0, 0, 0
    map_error, map_scaled_error, endpoint_error = 0., 0., 0.
    for path in FINAL.glob('*_models.npz'):
        selection = json.loads(path.with_name(path.name.replace('_models.npz', '_selection.json')).read_text())
        selected = selection['selected']; diagnostics = selection['diagnostics']
        for method, count in cfg['candidate_counts'].items():
            group = [c for c in selection['candidates'] if c['method'] == method]
            assert len(group) == count
            eligible = [c for c in group if c['valid']]
            minimum = min(c['score'] for c in eligible)
            eligible = [c for c in eligible if c['score']-minimum <= 1e-12*max(1., minimum)]
            expected = min(eligible, key=lambda c: (c['norm'], c['index']))
            assert selected[method]['index'] == expected['index']
            budget_checks += 1
        assert all(d['alpha_index'] == selected['S']['alpha_index'] for d in diagnostics)
        assert all(d['parameter'] == selected['S']['parameter'] for d in diagnostics)
        z = np.load(path); A, C, H = z['A'], z['Cemp'], z['CH1']
        R, centres, mc = z['R'], z['centres'], z['mc']
        sigma = float(path.stem.removesuffix('_models').rsplit('_', 1)[1]); variance = sigma*sigma/64
        gamma, alpha = selected['S']['gamma'], selected['S']['parameter']
        expected_maps = {
            0: plugin((1-gamma)*C+gamma*H, A, alpha, variance),
            1: plugin(C, A, selected['M0']['parameter'], variance),
            2: plugin(np.diag(np.diag(C)), A, selected['Diag']['parameter'], variance),
            5: plugin(H, A, selected['H1m']['parameter'], variance),
            6: plugin(C, A, alpha, variance),
            7: plugin(H, A, alpha, variance)}
        for j, expected in expected_maps.items():
            err = float(np.max(np.abs(expected-R[j])))
            scale = max(1., float(np.max(np.abs(R[j]))))
            map_error = max(map_error, err); map_scaled_error = max(map_scaled_error, err/scale)
            assert err/scale < 2e-7, (path.name, j, err/scale)
            np.testing.assert_allclose(centres[j], A@mc, rtol=1e-13, atol=1e-14)
        for j, alpha_endpoint in [(5, selected['H1m']['parameter']), (7, alpha)]:
            # Exact gamma=1 correspondence beta/s * Q, in coefficient space.
            h1 = np.linalg.solve(A.T@A+alpha_endpoint*variance/float(z['s'])*z['Q'], A.T)
            err = float(np.max(np.abs(h1-R[j]))); endpoint_error = max(endpoint_error, err)
            assert err/max(1., float(np.max(np.abs(R[j])))) < 2e-7
        if gamma == 0:
            np.testing.assert_array_equal(R[0], R[6]); np.testing.assert_array_equal(centres[0], centres[6])
            zero_count += 1
        for j in [3, 4]:
            np.testing.assert_array_equal(centres[j], z['my'])
        maps_count += 1

    # Direct integral and new seeds; 2 objects, 2 repetitions, 6 scenarios.
    scalar_checks, scalar_error, truth_error, empirical_map_error = 0, 0., 0., 0.
    for scenario, mask, p, height, partition in [
        ('primary', 'full', 2.5, .25, 30), ('partial', 'half', 2.5, .25, 30),
        ('rough_full', 'full', 1.5, .25, 31), ('rough_partial', 'half', 1.5, .25, 31),
        ('height_full', 'full', 2.5, .275, 30), ('height_partial', 'half', 2.5, .275, 30)]:
        z = np.load(FINAL / f'n64_g0_{mask}_0.01_models.npz')
        raw = np.load(FINAL / f'n64_g0_{scenario}_loss.npz')['loss']
        idx = np.arange(64 if mask == 'full' else 32)
        truth = np.load(FINAL / 'truth.npz')['rough' if p == 1.5 else 'nominal']
        for obj in [0, 399]:
            c = coefficients(partition, 0, obj, p)
            truth_error = max(truth_error, float(np.max(np.abs(c-truth[obj]))))
            np.testing.assert_allclose(c, truth[obj], rtol=1e-14, atol=1e-16)
            for rep in [0, 3]:
                noise = random(partition, 0, obj, 1, rep).normal(size=64)
                y = (analytic_transfer(height)@c+.01*noise)[idx]/8
                for j in range(8):
                    pred = z['mc']+z['R'][j]@(y-z['centres'][j])
                    error = np.r_[pred-c[:64], -c[64:]]
                    residual = z['A']@pred-y
                    values = np.array([error@error, error[8:]@error[8:], residual@residual])
                    scalar_error = max(scalar_error, float(np.max(np.abs(values-raw[obj, rep, j]))))
                    np.testing.assert_allclose(values, raw[obj, rep, j], rtol=2e-9, atol=1e-11)
                    scalar_checks += 3
        if scenario in ['primary', 'partial']:
            training = np.array([coefficients(10, 0, i) for i in range(64)])
            noise = np.array([random(10, 0, i, 1, 0).normal(size=64) for i in range(64)])
            y = (training@analytic_transfer(.25).T+.01*noise)[:, idx]/8
            c = training[:, :64]-training[:, :64].mean(axis=0); yy = y-y.mean(axis=0)
            cross, observed = c.T@yy/63, yy.T@yy/63
            selection = json.loads((FINAL / f'n64_g0_{mask}_0.01_selection.json').read_text())['selected']
            ridge = selection['E']['parameter'] * np.trace(observed)/len(idx)
            empirical = np.linalg.solve(observed+ridge*np.eye(len(idx)), cross.T).T
            empirical_map_error = max(empirical_map_error, float(np.max(np.abs(empirical-z['R'][3]))))
            np.testing.assert_allclose(empirical, z['R'][3], rtol=1e-7, atol=1e-8)

    # Descriptive only: min/median/max across the five training fits at N=64.
    secondary = []
    for scenario, *_ in cfg['scenarios']:
        entry = dict(scenario=scenario, N=64, summaries={})
        for denominator in ['M0', 'E', 'H1', 'Diag', 'M0_at_S', 'H1m_at_S']:
            ratios = np.array([float(table[64, g, scenario, 'S']['risk'])
                /float(table[64, g, scenario, denominator]['risk']) for g in range(5)])
            entry['summaries']['S/'+denominator] = dict(minimum=float(ratios.min()),
                median=float(np.median(ratios)), maximum=float(ratios.max()),
                ratios_above_one=int(np.count_nonzero(ratios > 1)))
        secondary.append(entry)
    report = dict(success=True, audit_type='independent archive recomputation; no campaign or selection rerun',
        frozen_hashes_verified=len(frozen['sha256']), source_hashes_verified=len(cfg['source_sha256'])+len(cfg['reused_source_sha256']),
        raw_archive_hashes_verified=len(analysis['raw_sha256']), summary_scalar_checks=checks,
        maximum_summary_difference=summary_error, failed_entries=failures,
        bootstrap_scalar_checks=int(saved_boot.size), maximum_bootstrap_difference=boot_error,
        primary=decisions, conclusion=conclusion,
        saved_model_files_checked=maps_count, selection_budget_and_rule_checks=budget_checks,
        gamma_zero_selected_identity_checks=zero_count,
        maximum_plugin_operator_difference=map_error,
        maximum_plugin_operator_scaled_difference=map_scaled_error,
        maximum_gamma_one_H1_endpoint_difference=endpoint_error,
        independently_reconstructed_empirical_maps=2, maximum_empirical_operator_difference=empirical_map_error,
        independent_scalar_loss_checks=scalar_checks, maximum_scalar_loss_difference=scalar_error,
        maximum_seeded_truth_difference=truth_error, secondary_descriptive_N64=secondary,
        interpretation='The useful 3% improvement is not established. The 3% partial-observation protection criterion is satisfied conditionally in each of five fits. This is not evidence of absence of every effect or general robustness.',
        limitations=['Approximate conditional percentile bootstrap, not finite-sample coverage guarantee.',
            'Five independent training fits share one physical test sample; no training-population confidence interval.',
            'S has 140 validation candidates versus 20 for M0; matched-alpha diagnostics separate covariance changes descriptively.',
            'Secondary comparisons are descriptive; nominal/model-error and roughness degradation must be retained.'])
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k != 'secondary_descriptive_N64'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', help='Fresh JSON destination; existing files are preserved')
    main(parser.parse_args().output)

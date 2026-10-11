"""Frozen-protocol analysis, integrity checked, without new scientific runs."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from laplace_learning.metrics import bootstrap_ratios, ordered_quantile, relative_errors
from laplace_learning.storage import METRICS, _atomic_bytes, atomic_npz, bounded_disk_usage, validate_condition


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1_000_000), b''):
            h.update(block)
    return h.hexdigest()


def clean(value):
    """Portable JSON: extended infinity is explicit, undefined remains null."""
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return '+inf' if value == np.inf else ('-inf' if value == -np.inf else 'undefined')
    return value


def summarize(values):
    return dict(median=ordered_quantile(values, .5), q95=ordered_quantile(values, .95),
                worst=ordered_quantile(values, 1.))


def select_diagnostics(errors, residuals, statuses):
    """Addendum 2: object means, deterministic ranks/ties, finite eligible case."""
    if errors.shape != residuals.shape or errors.ndim != 2:
        raise ValueError('object by repetition arrays required')
    e, r = np.mean(errors, axis=1), np.mean(residuals, axis=1)
    ids = np.arange(len(e))
    eligible = np.flatnonzero(np.isfinite(e) & np.isfinite(r) &
                              (r <= ordered_quantile(r, .1)) &
                              (e >= ordered_quantile(e, .9)))
    low = None if not len(eligible) else int(eligible[np.lexsort(
        (eligible, r[eligible], -e[eligible]))[0]])
    failed = np.argwhere(np.any(statuses != 0, axis=-1))
    return dict(median=int(np.flatnonzero(e == ordered_quantile(e, .5))[0]),
                q95=int(np.flatnonzero(e == ordered_quantile(e, .95))[0]),
                worst=int(np.flatnonzero(e == ordered_quantile(e, 1.))[0]),
                low_residual_high_error=low,
                failure=None if not len(failed) else failed[0].tolist(), repeat=0)


def run(source, output, prior_extra_cpu=0.):
    begin = time.process_time()
    source, output = Path(source), Path(output)
    report = json.loads((source/'report.json').read_text())
    if not report.get('complete') or report.get('verdict') != 'PASS':
        raise ValueError('full campaign report must be complete and PASS')
    if output.exists() and any(output.iterdir()):
        raise ValueError('preserve existing analysis; choose a new output directory')
    manifest = json.loads((source/'manifest.json').read_text())
    for name, meta in manifest.items():
        p = source/name
        if p.stat().st_size != meta['bytes'] or digest(p) != meta['sha256']:
            raise ValueError('campaign artifact integrity mismatch: '+name)
    config = json.loads((source/'configuration.json').read_text())
    cfg = config['design']
    state = json.loads((source/'state.json').read_text())
    arrays = {k: np.load(source/(k+'.npy'), mmap_mode='r', allow_pickle=False)
              for k in METRICS+('status',)}
    expected = {f'{g}_{m}_{n}' for g in range(cfg['games'])
                for m in range(len(cfg['masks'])) for n in range(len(cfg['levels']))}
    if set(state['completed']) != expected or set(state['condition_hashes']) != expected:
        raise ValueError('missing completed conditions')
    for key in sorted(expected):
        validate_condition(arrays, tuple(map(int, key.split('_'))), state['condition_hashes'][key])
    truth_meta = json.loads((source/'truth_manifest.json').read_text())
    truth = {}
    for name, chunks in truth_meta.items():
        parts = []
        for chunk in chunks:
            p = source/chunk['path']
            if digest(p) != chunk['sha256']:
                raise ValueError('truth chunk hash mismatch')
            parts.append(np.load(p, allow_pickle=False))
        truth[name] = np.concatenate(parts)
    norms = {name: np.linalg.norm(a, axis=1) for name, a in truth.items()}
    methods = config['methods']
    records, primary, draws, diagnostics = [], [], {}, {}
    for g in range(cfg['games']):
        d0 = 1e-12*np.median(norms[f'train_{g}'])
        for mi, mask in enumerate(cfg['masks']):
            for ni, nu in enumerate(cfg['levels']):
                selected = json.loads((source/f'selection_{g}_{mi}_{ni}.json').read_text())
                candidates = json.loads((source/f'candidates_{g}_{mi}_{ni}.json').read_text())
                for ri, regime in enumerate(cfg['regimes']):
                    errors = relative_errors(arrays['error_abs'][g, mi, ni, ri],
                        norms['rough' if regime == 'rough' else 'nominal'][:, None, None], d0)
                    statuses = arrays['status'][g, mi, ni, ri]
                    for j, method in enumerate(methods):
                        status = statuses[..., j]
                        grid = [c['grid_index'] for c in candidates if c['method'] == method]
                        on_edge = len(grid) > 1 and selected[method]['grid_index'] in (min(grid), max(grid))
                        row = dict(training=g, mask=mask, noise=nu, regime=regime, method=method,
                                   count=int(status.size), failures=int(np.count_nonzero(status)),
                                   numerical=int(np.count_nonzero(status & 1)),
                                   metric=int(np.count_nonzero(status & 2)),
                                   family=int(np.count_nonzero(status & 4)),
                                   frequency_gt_one=float(np.mean(errors[..., j] > 1)),
                                   operator_norm=selected[method]['operator_norm'],
                                   grid_index=selected[method]['grid_index'],
                                   parameter=selected[method]['parameter'], grid_edge=on_edge)
                        for label, data in [('relative', errors[..., j])] + [
                                (k, arrays[k][g, mi, ni, ri, ..., j]) for k in METRICS]:
                            row.update({label+'_'+stat: v for stat, v in summarize(data).items()})
                        records.append(row)
        for contrast, mask, regime in [('P1', 'contiguous', 'nominal'), ('P2', 'full', 'height_5')]:
            mi, ni, ri = cfg['masks'].index(mask), cfg['levels'].index(.01), cfg['regimes'].index(regime)
            rn = cfg['regimes'].index('nominal')
            cap, base = methods.index('learned_signed_cap'), methods.index('learned_signed')
            def err(r, j):
                return relative_errors(arrays['error_abs'][g, mi, ni, r, ..., j], norms['nominal'], d0)
            result = bootstrap_ratios(err(ri, cap), err(ri, base), err(rn, cap), err(rn, base), g)
            draws[f'{contrast}_training_{g}'] = result.pop('draws')
            primary.append(dict(contrast=contrast, training=g, **result))
            if g == 0:
                diagnostics[contrast] = select_diagnostics(err(ri, base),
                    arrays['residual_nominal'][g, mi, ni, ri, ..., base],
                    arrays['status'][g, mi, ni, ri])
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
    csv_bytes = stream.getvalue().encode('utf-8')
    primary_data = dict(comparisons=primary, diagnostics=diagnostics,
        robust_across_trainings=all(p['verdict'] == 'SATISFIED' for p in primary),
        coverage='Bonferroni four ratios conditional on each training; no simultaneous three-training coverage',
        quantile='ceil(q*N), no interpolation', bootstrap_object_draws=2000,
        campaign_manifest_sha256=digest(source/'manifest.json'))
    primary_bytes = (json.dumps(clean(primary_data), ensure_ascii=False, indent=2, allow_nan=False)+'\n').encode()
    buffer = io.BytesIO()
    np.savez(buffer, **draws)
    estimated = len(csv_bytes)+len(primary_bytes)+len(buffer.getvalue())+10_000
    if estimated > 1_000_000:
        raise RuntimeError('analysis group exceeds 1 MB allocation')
    bounded_disk_usage(source.parent, estimated+1_000_000)
    output.mkdir(parents=True, exist_ok=True)
    _atomic_bytes(output/'summary.csv', csv_bytes)
    _atomic_bytes(output/'primary.json', primary_bytes)
    atomic_npz(output/'bootstrap.npz', **draws)
    used = time.process_time()-begin
    cumulative = report['cpu_seconds']+prior_extra_cpu+used
    if cumulative > 14400 or used > 1800:
        raise RuntimeError('analysis CPU allocation exceeded')
    evidence = dict(verdict='PASS', complete=True, campaign_integrity_verified=True, conditions=len(expected),
        aggregate_rows=len(records), failures=sum(r['failures'] for r in records),
        analysis_cpu_seconds=used, cpu_seconds=cumulative, cumulative_cpu_seconds=cumulative, extra_prior_cpu_seconds=prior_extra_cpu,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        files={p.name: dict(bytes=p.stat().st_size, sha256=digest(p)) for p in output.iterdir()},
        analysis_script_sha256=digest(__file__))
    _atomic_bytes(output/'report.json', (json.dumps(evidence, indent=2)+'\n').encode())
    if bounded_disk_usage(output) > 1_000_000:
        raise RuntimeError('analysis group byte allocation exceeded')
    return evidence, primary_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', default='results/full')
    parser.add_argument('--output', default='results/analysis')
    parser.add_argument('--prior-extra-cpu', type=float, default=0.)
    args = parser.parse_args()
    report, primary = run(args.source, args.output, args.prior_extra_cpu)
    print(json.dumps(clean(dict(report=report, primary=primary)), ensure_ascii=False, indent=2))

"""Frozen affine estimators and observation-only model selection (protocol §§3–8).

Samples are rows. A maps physical coefficients to normalized observations.
Invalid candidates are retained with explicit failure metadata; no method fallback.
"""
import numpy as np
from scipy.linalg import cho_factor, cho_solve

METHODS = (
    'tikhonov_l2', 'tikhonov_h1', 'tsvd', 'wiener', 'wiener_cap',
    'learned_signed', 'learned_signed_cap', 'learned_positive',
    'learned_bounded', 'learned_bounded_cap', 'lmmse',
)


def _solve_spd(matrix, rhs):
    return cho_solve(cho_factor(matrix, lower=True, check_finite=True), rhs)


def fit_candidates(A, c_train, d_train, noise_variance):
    """Return <=144 candidate dictionaries, fitting only supplied training pairs.

    Nonfinite input is a data error (ValueError). Numerical candidate failures
    remain in the returned list with valid=False and a failure description.
    """
    A, c, d = (np.asarray(x, dtype=float) for x in (A, c_train, d_train))
    if A.ndim != 2 or c.ndim != 2 or d.ndim != 2:
        raise ValueError('A, c_train and d_train must be matrices')
    m, n = A.shape
    if not m or not n or c.shape[0] < 2 or c.shape != (d.shape[0], n) or d.shape[1] != m:
        raise ValueError('incompatible dimensions or fewer than two training objects')
    if not all(np.isfinite(x).all() for x in (A, c, d)):
        raise ValueError('nonfinite training data or operator')
    if not np.isfinite(noise_variance) or noise_variance < 0:
        raise ValueError('noise_variance must be finite and nonnegative')
    mu_c, mu_d = c.mean(axis=0), d.mean(axis=0)
    cc, dd = c - mu_c, d - mu_d
    divisor = c.shape[0] - 1
    u, singular, vh = np.linalg.svd(A, full_matrices=False)
    rank = int(np.sum(singular > 1e-12 * singular[0]))
    ur, vr, sr = u[:, :rank], vh[:rank].T, singular[:rank]
    beta, z = cc @ vr, dd @ ur
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        vz = np.sum(z * z, axis=0) / divisor
        p = np.sum(beta * beta, axis=0) / divisor
        covariance = np.sum(beta * z, axis=0) / divisor
        threshold = 100 * np.finfo(float).eps * max(1., float(np.max(vz, initial=0.)))
        learned = np.zeros(rank)
        retained = vz > threshold
        # A nonfinite variance is a numerical failure, never a zero-gain rejection.
        retained |= ~np.isfinite(vz)
        learned[retained] = covariance[retained] / vz[retained]
        denominator = sr * sr * p + noise_variance
        wiener = np.divide(sr * p, denominator, out=np.zeros(rank), where=denominator != 0)
        Ccd, Cdd = cc.T @ dd / divisor, dd.T @ dd / divisor
    # Preserve numerical failures through clipping: an infinite base gain must
    # not become a valid finite capped estimator. NaN marks an invalid family.
    bad_learned = (~np.isfinite(learned) | ~np.isfinite(covariance) | ~np.isfinite(vz))
    learned[bad_learned] = np.nan
    wiener[~np.isfinite(wiener) | ~np.isfinite(p) | ~np.isfinite(denominator)] = np.nan
    scale = max(float(np.trace(Cdd) / m), 1e-15)
    candidates = []

    def add(method, parameter, grid_index, make):
        record = dict(method=method, parameter=parameter, grid_index=int(grid_index),
                      mu_c=mu_c.copy(), mu_d=mu_d.copy(), numerical_rank=rank,
                      variance_threshold=float(threshold),
                      rejected_variances=int(np.sum(vz <= threshold)))
        try:
            with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
                R = make()
            if not np.isfinite(R).all():
                raise FloatingPointError('nonfinite reconstruction matrix')
            norm = float(np.linalg.norm(R, ord=2))
            if not np.isfinite(norm):
                raise FloatingPointError('nonfinite operator norm')
            record.update(R=R, operator_norm=norm, valid=True, failure=None)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as exc:
            record.update(R=np.full((n, m), np.nan), operator_norm=np.inf,
                          valid=False, failure=f'{type(exc).__name__}: {exc}')
        candidates.append(record)

    def spectral(g):
        return (vr * g) @ ur.T

    alphas, caps = np.logspace(-10, 1, 20), np.logspace(0, 4, 20)
    penalty = 1 + (np.arange(n) * np.pi) ** 2
    penalty /= penalty.max()
    for index, alpha in enumerate(alphas):
        # Full normal equations, including numerical SVD-discarded directions.
        add('tikhonov_l2', float(alpha), index,
            lambda: _solve_spd(A.T @ A + alpha * np.eye(n), A.T))
        add('tikhonov_h1', float(alpha), index,
            lambda: _solve_spd(A.T @ A + alpha * np.diag(penalty), A.T))
        add('lmmse', float(alpha), index,
            lambda: _solve_spd(Cdd + alpha * scale * np.eye(m), Ccd.T).T)
    ranks = np.unique(np.rint(np.linspace(1, rank, 20)).astype(int)) if rank else [0]
    for index, k in enumerate(ranks):
        add('tsvd', int(k), index, lambda: spectral(np.where(np.arange(rank) < k, 1 / sr, 0)))
    fixed = [('wiener', wiener), ('learned_signed', learned),
             ('learned_positive', np.maximum(0, learned)),
             ('learned_bounded', np.clip(learned, 0, 1 / sr))]
    for method, gains in fixed:
        add(method, None, 0, lambda: spectral(gains))
    for index, cap in enumerate(caps):
        add('wiener_cap', float(cap), index, lambda: spectral(np.minimum(wiener, cap)))
        add('learned_signed_cap', float(cap), index, lambda: spectral(np.clip(learned, -cap, cap)))
        add('learned_bounded_cap', float(cap), index,
            lambda: spectral(np.clip(learned, 0, np.minimum(1 / sr, cap))))
    return candidates


def predict(estimator, d):
    """Apply affine estimator to a row batch (or one observation vector).

    An invalid estimator produces NaNs, making its failure explicit downstream.
    """
    observations = np.asarray(d, dtype=float)
    with np.errstate(over='ignore', invalid='ignore'):
        return estimator['mu_c'] + (observations - estimator['mu_d']) @ estimator['R'].T


def select_candidates(candidates, A, dval1, dval2, *, return_diagnostics=False):
    """Select by independent-repetition predictive risk, without validation truth.

    Invalid predictions/scores get +inf. Ties within the protocol tolerance are
    resolved by the true matrix operator norm and then fixed grid index.
    return_diagnostics=True also returns scalar metadata for every candidate,
    without retaining reconstruction matrices or reevaluating validation.
    """
    A, first, second = (np.asarray(x, dtype=float) for x in (A, dval1, dval2))
    if first.ndim != 2 or first.shape != second.shape or first.shape[0] == 0 or first.shape[1] != A.shape[0]:
        raise ValueError('validation batches must be nonempty compatible row matrices')
    groups = {}
    diagnostics = []
    for candidate in candidates:
        result = dict(candidate)
        score = np.inf
        if candidate['valid']:
            reconstruction = predict(candidate, first)
            with np.errstate(over='ignore', invalid='ignore'):
                residual = reconstruction @ A.T - second
                value = np.mean(np.sum(residual * residual, axis=1))
            if np.isfinite(reconstruction).all() and np.isfinite(value):
                score = float(value)
        result['validation_score'] = score
        result['score'] = score
        result['valid'] = bool(candidate['valid'] and np.isfinite(score))
        if candidate['valid'] and not result['valid']:
            result['failure'] = 'nonfinite_validation_prediction_or_score'
        diagnostics.append({key: result[key] for key in
                            ('method', 'parameter', 'grid_index', 'operator_norm',
                             'score', 'valid', 'failure') if key in result})
        groups.setdefault(candidate['method'], []).append(result)
    selected = {}
    for method, group in groups.items():
        valid = [x for x in group if x['valid']]
        if not valid:
            result = dict(group[0])
            result.update(valid=False, failure='all_candidates_invalid',
                          parameter=None, grid_index=-1,
                          R=np.full_like(result['R'], np.nan), operator_norm=np.inf)
        else:
            minimum = min(x['score'] for x in valid)
            tolerance = 1e-12 * max(1., abs(minimum))
            tied = [x for x in valid if x['score'] - minimum <= tolerance]
            result = min(tied, key=lambda x: (x['operator_norm'], x['grid_index']))
        selected[method] = result
    return (selected, diagnostics) if return_diagnostics else selected

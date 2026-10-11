"""Frozen study 3 covariance families; no access to validation or test in fitting."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve


METHODS = ['S', 'M0', 'Diag', 'E', 'H1', 'H1m']
FIXED = ['M0_at_S', 'H1m_at_S']
GAMMAS = np.array([0., .01, .03, .1, .3, .6, 1.])
ALPHAS = np.logspace(-3, 3, 20)
EMPIRICAL_LAMBDAS = np.logspace(-8, 1, 20)
H1_PARAMETERS = np.logspace(-10, 1, 20)


def _solve(matrix, rhs):
    """Solve the stated SPD system, without inverse, jitter, or fallback."""
    return cho_solve(cho_factor(matrix, lower=True), rhs)


def fit_candidates(A, c_train, y_train, noise_variance):
    """Fit the protocol's 240 candidates from training arrays only.

    ``noise_variance`` is sigma**2 / 64 after fixed observation scaling.
    Q is returned as a matrix, in the unnormalised H1 convention.
    A candidate's ``parameter`` is alpha, lambda, or historical H1 a.
    ``alpha_index`` is meaningful for S/M0/Diag/H1m only (None otherwise).
    """
    A = np.asarray(A, dtype=float)
    c = np.asarray(c_train, dtype=float)
    y = np.asarray(y_train, dtype=float)
    if (A.ndim != 2 or c.ndim != 2 or y.ndim != 2 or min(A.shape) < 1
            or c.shape[0] < 2 or c.shape[0] != y.shape[0]
            or c.shape[1] != A.shape[1] or y.shape[1] != A.shape[0]):
        raise ValueError('incompatible training shapes')
    if not all(np.isfinite(x).all() for x in (A, c, y)):
        raise ValueError('nonfinite training input')
    if np.ndim(noise_variance) != 0:
        raise ValueError('noise variance must be a positive finite scalar')
    noise_variance = float(noise_variance)
    if not np.isfinite(noise_variance) or noise_variance <= 0:
        raise ValueError('noise variance must be a positive finite scalar')

    mc, my = c.mean(0), y.mean(0)
    cc, yy = c - mc, y - my
    den = c.shape[0] - 1
    Cemp = cc.T @ cc / den
    Ccy = cc.T @ yy / den
    Cyy = yy.T @ yy / den
    energy = float(np.trace(Cemp))
    ty = float(np.trace(Cyy) / A.shape[0])
    if not np.isfinite(energy) or energy <= 0:
        raise ValueError('signal covariance trace must be positive and finite')
    if not np.isfinite(ty) or ty <= 0:
        raise ValueError('observation covariance scale must be positive and finite')
    q = 1 + (np.pi * np.arange(A.shape[1], dtype=float))**2
    Q = np.diag(q)
    s = float(energy / np.sum(1 / q))
    CH1 = np.diag(s / q)
    model_centre = A @ mc
    eye = np.eye(A.shape[0])
    candidates = []

    def add(method, index, alpha_index, parameter, gamma, centre, make):
        try:
            R = make()
            if not np.isfinite(R).all():
                raise FloatingPointError('nonfinite reconstruction operator')
            norm = float(np.linalg.norm(R, 2))
            if not np.isfinite(norm):
                raise FloatingPointError('nonfinite operator norm')
            valid, failure = True, None
        except (np.linalg.LinAlgError, ValueError, FloatingPointError) as exc:
            R = np.full((A.shape[1], A.shape[0]), np.nan)
            norm, valid, failure = float('inf'), False, str(exc)
        candidates.append(dict(
            method=method, index=int(index), alpha_index=alpha_index,
            parameter=float(parameter), gamma=gamma, R=R, centre=centre.copy(),
            norm=norm, valid=valid, failure=failure))

    def plugin(C, alpha):
        CA = C @ A.T
        return _solve(A @ CA + alpha * noise_variance * eye, CA.T).T

    for gamma_index, gamma in enumerate(GAMMAS):
        Cgamma = (1 - gamma) * Cemp + gamma * CH1
        for alpha_index, alpha in enumerate(ALPHAS):
            add('S', gamma_index * len(ALPHAS) + alpha_index, alpha_index,
                alpha, float(gamma), model_centre,
                lambda: plugin(Cgamma, alpha))
    for method, C, gamma in [('M0', Cemp, 0.),
                             ('Diag', np.diag(np.diag(Cemp)), None),
                             ('H1m', CH1, 1.)]:
        for alpha_index, alpha in enumerate(ALPHAS):
            add(method, alpha_index, alpha_index, alpha, gamma, model_centre,
                lambda: plugin(C, alpha))
    for index, parameter in enumerate(EMPIRICAL_LAMBDAS):
        add('E', index, None, parameter, None, my,
            lambda: _solve(Cyy + parameter * ty * eye, Ccy.T).T)
    normal = A.T @ A
    Qn = np.diag(q / q.max())
    for index, parameter in enumerate(H1_PARAMETERS):
        add('H1', index, None, parameter, None, my,
            lambda: _solve(normal + parameter * Qn, A.T))
    return dict(candidates=candidates, mc=mc, my=my, Cemp=Cemp,
                CH1=CH1, s=s, Q=Q)


def predict(family, candidate, y):
    """Apply the candidate's own centre (A*mc or empirical my)."""
    y = np.asarray(y, dtype=float)
    return family['mc'] + (y - candidate['centre']) @ candidate['R'].T


def select(family, yval, cval128):
    """Select on explicitly supplied validation data, preserving all failures."""
    y = np.asarray(yval, dtype=float)
    truth = np.asarray(cval128, dtype=float)
    n = len(family['mc'])
    m = len(family['my'])
    if (y.ndim != 2 or truth.ndim != 2 or len(y) == 0
            or len(y) != len(truth) or y.shape[1] != m or truth.shape[1] < n):
        raise ValueError('incompatible validation shapes')
    if not np.isfinite(y).all() or not np.isfinite(truth).all():
        raise ValueError('nonfinite validation input')
    omitted = np.sum(truth[:, n:]**2, axis=1)
    scores = []
    for source in family['candidates']:
        candidate = dict(source)
        score = float('inf')
        if candidate['valid']:
            prediction = predict(family, candidate, y)
            score = float(np.mean(np.sum((prediction - truth[:, :n])**2, axis=1)
                                  + omitted))
        candidate['score'] = score
        if not np.isfinite(score):
            candidate.update(valid=False,
                             failure=candidate['failure'] or 'nonfinite validation loss')
        scores.append(candidate)
    selected = {}
    for method in METHODS:
        all_candidates = [candidate for candidate in scores if candidate['method'] == method]
        if not all_candidates:
            raise ValueError('missing candidate family: ' + method)
        valid = [candidate for candidate in all_candidates if candidate['valid']]
        if not valid:
            # A record of a failed selection, never a different method or operator.
            selected[method] = dict(all_candidates[0], valid=False, score=float('inf'))
            continue
        minimum = min(candidate['score'] for candidate in valid)
        tolerance = 1e-12 * max(1., minimum)
        selected[method] = min(
            (candidate for candidate in valid if candidate['score'] - minimum <= tolerance),
            key=lambda candidate: (candidate['norm'], candidate['index']))
    fixed = []
    for method, diagnostic in zip(['M0', 'H1m'], FIXED):
        candidate = dict(next(candidate for candidate in scores
                              if candidate['method'] == method
                              and candidate['alpha_index'] == selected['S']['alpha_index']))
        candidate['diagnostic'] = diagnostic
        if not selected['S']['valid']:
            candidate.update(valid=False, score=float('inf'),
                             R=np.full_like(candidate['R'], np.nan),
                             failure='S selection invalid; matched-alpha diagnostic unavailable')
        fixed.append(candidate)
    return selected, fixed, scores

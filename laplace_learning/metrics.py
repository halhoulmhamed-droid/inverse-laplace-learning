"""Protocol A metrics: extended-real errors and paired object bootstrap.

Missing/corrupt input is rejected, not confused with a failed reconstruction.
Status bits: 1 = nonfinite reconstruction; 2 = metric overflow.
"""
import math

import numpy as np


def _errors(values):
    if values is None:
        raise ValueError("missing errors")
    a = np.asarray(values, dtype=float)
    if not a.size or np.any(np.isnan(a)) or np.any(a < 0):
        raise ValueError("errors must be nonempty, nonnegative and without NaN")
    return a


def ordered_quantile(values, q):
    """Nearest rank ceil(q*N), with +inf retained and no interpolation."""
    if not 0 < q <= 1:
        raise ValueError("q must be in (0, 1]")
    a = _errors(values).ravel()
    k = math.ceil(q * a.size) - 1
    return float(np.partition(a, k)[k])


def relative_errors(absolute, true_norm, d0):
    a = _errors(absolute)
    t = np.asarray(true_norm, dtype=float)
    if not np.isfinite(d0) or d0 <= 0 or np.any(~np.isfinite(t)) or np.any(t < 0):
        raise ValueError("finite nonnegative truth norms and positive d0 required")
    if a.ndim == 2 and t.shape == (a.shape[0],):
        t = t[:, None]
    with np.errstate(over="ignore"):
        return a / np.maximum(t, d0)


def paired_ratio_summary(num, den):
    """Ratio of supplied scalar summaries; no victory over an infinite baseline."""
    num, den = float(num), float(den)
    if np.isnan(num) or num < 0 or not np.isfinite(den) or den <= 0:
        return {"ratio": None, "defined": False, "verdict": "NON_CONCLUSIVE"}
    with np.errstate(over="ignore"):
        ratio = float(np.divide(num, den))
    return {"ratio": ratio, "defined": True, "verdict": "DEFINED"}


def bootstrap_ratios(errors_capped, errors_base, nominal_capped, nominal_base,
                     train_id, n_boot=2000):
    """Conditional paired bootstrap; rows are objects, columns repetitions.

    Call separately for the two contrasts with the same train_id to reuse
    exactly the same object draws for all four Bonferroni-corrected ratios.
    An undefined draw invalidates its bound; it is never removed.
    """
    arrays = [_errors(x) for x in
              (errors_capped, errors_base, nominal_capped, nominal_base)]
    if any(x.ndim != 2 or x.shape != arrays[0].shape for x in arrays):
        raise ValueError("four matching object-by-repetition arrays required")
    if int(train_id) != train_id or train_id < 0 or int(n_boot) != n_boot or n_boot < 1:
        raise ValueError("nonnegative integer train_id and positive n_boot required")

    def ratios(xs):
        return [paired_ratio_summary(ordered_quantile(xs[0], .95),
                                     ordered_quantile(xs[1], .95))["ratio"],
                paired_ratio_summary(ordered_quantile(xs[2], .5),
                                     ordered_quantile(xs[3], .5))["ratio"]]

    point = ratios(arrays)
    draws = np.full((n_boot, 2), np.nan)
    for b in range(n_boot):
        rng = np.random.Generator(np.random.PCG64(
            np.random.SeedSequence([20261011, 5, int(train_id), 0, 2, b])))
        ids = rng.integers(0, arrays[0].shape[0], size=arrays[0].shape[0])
        draw = ratios([x[ids] for x in arrays])
        draws[b] = [np.nan if x is None else x for x in draw]
    invalid = np.isnan(draws).sum(axis=0)
    upper = [None if invalid[j] else ordered_quantile(draws[:, j], .9875)
             for j in range(2)]
    if any(x is None for x in point + upper):
        verdict = "NON_CONCLUSIVE"
    elif all(np.isfinite(x) for x in point + upper) and upper[0] < .90 and upper[1] < 1.05:
        verdict = "SATISFIED"
    else:
        verdict = "NOT_SATISFIED"
    return dict(q95_ratio=point[0], median_ratio=point[1],
                upper_q95=upper[0], upper_median=upper[1],
                invalid_q95=int(invalid[0]), invalid_median=int(invalid[1]),
                draws=draws, verdict=verdict)


def _row_norm(a):
    """Scaled Euclidean norm avoids squaring large finite coefficients."""
    if a.shape[1] == 0:
        return np.zeros(a.shape[0])
    scale = np.max(np.abs(a), axis=1)
    safe = np.where((scale > 0) & np.isfinite(scale), scale, 1.)
    with np.errstate(over="ignore", invalid="ignore"):
        return safe * np.sqrt(np.sum((a / safe[:, None]) ** 2, axis=1))


def metrics_for_predictions(chat, ctrue_full, A, d, repeated_d, Vr):
    """Five rowwise metrics, including tail only in total error.

    Vr contains orthonormal retained right singular vectors as columns.
    A has shape (observations, inverse coefficients). All other arrays use
    samples as rows. Truth and observation corruption is an input error.
    """
    args = (chat, ctrue_full, A, d, repeated_d, Vr)
    if any(x is None for x in args):
        raise ValueError("missing input: campaign is incomplete")
    chat, truth, A, d, repeated_d, Vr = [np.asarray(x, dtype=float) for x in args]
    if any(x.ndim != 2 for x in (chat, truth, A, d, repeated_d, Vr)):
        raise ValueError("all arrays must be two-dimensional")
    rows, n = chat.shape
    if (truth.shape[0] != rows or truth.shape[1] < n or A.shape[1] != n
            or d.shape != (rows, A.shape[0]) or repeated_d.shape != d.shape
            or Vr.shape[0] != n):
        raise ValueError("incompatible metric dimensions")
    if any(np.any(~np.isfinite(x)) for x in (truth, A, d, repeated_d, Vr)):
        raise ValueError("nonfinite source input: corrupt data, not method failure")
    numerical = np.any(~np.isfinite(chat), axis=1)
    safe_chat = np.where(numerical[:, None], 0., chat)
    with np.errstate(over="ignore", invalid="ignore"):
        delta = safe_chat - truth[:, :n]
        visible = delta @ Vr
        null = delta - visible @ Vr.T
        tail = truth[:, n:]
        prediction = safe_chat @ A.T
        values = np.column_stack([
            _row_norm(np.concatenate((delta, tail), axis=1)),
            _row_norm(visible),
            _row_norm(null),
            _row_norm(prediction - d),
            _row_norm(prediction - repeated_d)])
    metric_failure = np.any(~np.isfinite(values), axis=1) & ~numerical
    status = numerical.astype(np.uint8) | (metric_failure.astype(np.uint8) << 1)
    values[status != 0] = np.inf
    names = ("error_abs", "visible", "null", "residual_nominal", "residual_independent")
    return {**{name: values[:, j] for j, name in enumerate(names)}, "status": status}

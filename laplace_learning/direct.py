"""Independent harmonic and reflected finite-difference forward operators.

Columns are orthonormal cosine coefficients; rows are physical interval means.
No mask, weighting or normalization is applied implicitly. Grids use numbers of
intervals, not numbers of points. The numerical solver really factors a sparse
2-D mixed-boundary system; modal formulas are only independent diagnostics.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


def _positive_int(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _height(H):
    if not np.isfinite(H) or H <= 0:
        raise ValueError("H must be finite and positive")
    return float(H)


def cosine_sensor_matrix(K, sensors=64):
    """Return exact interval means of phi_0=1, phi_k=sqrt(2)cos(k*pi*x).

    The equivalent midpoint/sinc identity avoids subtraction of nearby sines.
    Each sensor has width 1/sensors; this is a mean, not an integral.
    """
    K, sensors = _positive_int(K, "K"), _positive_int(sensors, "sensors")
    k = np.arange(K, dtype=float)
    midpoints = (np.arange(sensors) + 0.5) / sensors
    matrix = np.sqrt(2.0) * np.cos(np.pi * midpoints[:, None] * k)
    matrix *= np.sinc(k / (2 * sensors))
    matrix[:, 0] = 1.0
    return matrix


def analytic_matrix(K, H=0.25, sensors=64):
    """Exact harmonic transfer followed by exact physical sensor averages."""
    H = _height(H)
    matrix = cosine_sensor_matrix(K, sensors)
    x = np.pi * H * np.arange(matrix.shape[1])
    decay = np.exp(-x)
    return matrix * (2 * decay / (1 + decay * decay))


def mask_indices(name, sensors=64):
    """Known masks; observed weights remain 1/sensors after selection."""
    sensors = _positive_int(sensors, "sensors")
    if name == "full":
        return np.arange(sensors)
    if name == "even":
        return np.arange(0, sensors, 2)
    if name == "contiguous":
        if sensors % 2:
            raise ValueError("contiguous half-mask requires an even sensor count")
        return np.arange(sensors // 2)
    raise ValueError("mask must be full, even or contiguous")


def piecewise_linear_sensor_matrix(nx, sensors=64):
    """Exactly integrate the nodal linear interpolant, including sub-cell sensors.

    Returns sensors by (nx+1); multiplication by nodal values gives means.
    This acquisition matrix alone has rank 17 on the first 32 sensors when
    nx=32 and sensors=64. It says nothing by itself about numerical PDE rank.
    """
    nx, sensors = _positive_int(nx, "nx"), _positive_int(sensors, "sensors")
    weights = np.zeros((sensors, nx + 1))
    for j in range(sensors):
        a, b = j / sensors, (j + 1) / sensors
        for i in range(max(0, int(np.floor(a * nx))), min(nx, int(np.ceil(b * nx)))):
            left, right = max(a, i / nx), min(b, (i + 1) / nx)
            if right <= left:
                continue
            # Integrals of local hat functions, scaled to sensor averages.
            width = (right - left) * sensors
            local_midpoint = nx * ((left + right) / 2) - i
            weights[j, i] += width * (1 - local_midpoint)
            weights[j, i + 1] += width * local_midpoint
    return weights


def _reflected_1d(intervals, spacing):
    """Centered second derivative including reflected boundary ghost nodes."""
    length = intervals + 1
    derivative = sparse.diags([np.ones(length - 1), -2 * np.ones(length),
                               np.ones(length - 1)], [-1, 0, 1], format="lil")
    derivative[0, 1] = 2.0
    derivative[-1, -2] = 2.0
    return derivative.tocsc() / spacing**2


def reflected_laplacian(nx, ny, H=0.25):
    """Return (B, D) for B*u=D*q with top Dirichlet and reflected other edges.

    Unknowns are ordered y-major, x-minor, including y=0 but excluding y=H.
    B is minus the discrete Laplacian; D injects the top nodal values. The
    reflection stencil is not symmetric in unweighted Euclidean coordinates.
    """
    nx, ny, H = _positive_int(nx, "nx"), _positive_int(ny, "ny"), _height(H)
    lx = _reflected_1d(nx, 1 / nx)
    ly = _reflected_1d(ny, H / ny)
    B = -(sparse.kron(sparse.eye(ny), lx, format="csc") +
          sparse.kron(ly[:ny, :ny], sparse.eye(nx + 1), format="csc"))
    D = sparse.kron(ly[:ny, ny:], sparse.eye(nx + 1), format="csc")
    return B, D


def fd_matrix(n, nx=64, ny=16, H=0.25, sensors=64):
    """Solve the sparse two-dimensional PDE and average its lower trace.

    n <= nx+1 prevents aliased top-boundary modes. Eight right hand sides at
    a time bound work memory; interior fields are discarded after acquisition.
    """
    n, nx, ny = _positive_int(n, "n"), _positive_int(nx, "nx"), _positive_int(ny, "ny")
    sensors, H = _positive_int(sensors, "sensors"), _height(H)
    if n > nx + 1:
        raise ValueError("n cannot exceed the number of horizontal nodes")
    # Guard sparse factorization against accidental grids outside this study.
    if (nx + 1) * ny > 100_000:
        raise ValueError("grid exceeds the supported bounded-memory study size")
    B, D = reflected_laplacian(nx, ny, H)
    factor = splu(B)
    acquisition = piecewise_linear_sensor_matrix(nx, sensors)
    result = np.empty((sensors, n))
    x = np.linspace(0, 1, nx + 1)
    for start in range(0, n, 8):
        modes = np.arange(start, min(start + 8, n))
        top = np.sqrt(2.0) * np.cos(np.pi * x[:, None] * modes)
        if start == 0:
            top[:, 0] = 1.0
        field = factor.solve(np.asarray(D @ top))
        result[:, start:start + len(modes)] = acquisition @ field[:nx + 1]
    if not np.isfinite(result).all():
        raise FloatingPointError("non-finite finite-difference solution")
    return result


def modal_discrete_ratio(k, nx, ny, H=0.25):
    """Discrete/continuous sensor-vector ratio for grid-aligned sensor edges.

    Only valid for 0 <= k < nx and sensors dividing nx; the sensor count does
    not enter the resulting scalar formula. Do not use on the nx=32 inverse
    with 64 sensors. This diagnostic is never used by the sparse PDE solver.
    """
    nx, ny, H = _positive_int(nx, "nx"), _positive_int(ny, "ny"), _height(H)
    if not isinstance(k, (int, np.integer)) or not 0 <= k < nx:
        raise ValueError("k must be an integer in [0,nx)")
    if k == 0:
        return 1.0
    s = k * np.pi / nx
    gamma = 2 * np.arcsinh((H * nx / ny) * np.sin(s / 2))
    # log(cosh) avoids overflow before cancellation of close large exponents.
    logcosh = lambda z: np.logaddexp(z, -z) - np.log(2.0)
    return float(np.exp(logcosh(k * np.pi * H) - logcosh(ny * gamma)) *
                 ((s / 2) / np.tan(s / 2)))

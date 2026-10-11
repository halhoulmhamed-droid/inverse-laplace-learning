"""Bounded atomic artifacts and integrity checks for completed conditions."""
import hashlib
import io
import os
from pathlib import Path
import tempfile
import numpy as np

MAX_TRANSIENT = 1_000_000
METRICS = ('error_abs', 'visible', 'null', 'residual_nominal', 'residual_independent')


def bounded_disk_usage(root, planned_bytes=0, limit=200_000_000):
    """Count all files below the shared experiment root, including temporaries.

    planned_bytes must include headers and any reserved yet unallocated budget.
    The caller supplies the common root covering pilots, campaigns and analysis.
    """
    if int(planned_bytes) != planned_bytes or planned_bytes < 0 or limit < 0:
        raise ValueError('invalid byte budget')
    used = sum(p.stat().st_size for p in Path(root).rglob('*') if p.is_file())
    if used + planned_bytes > limit:
        raise RuntimeError('planned storage exceeds the shared output budget')
    return used


def _atomic_bytes(path, data, max_bytes=MAX_TRANSIENT):
    if len(data) > max_bytes:
        raise ValueError('atomic temporary exceeds its byte limit; split the array')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def atomic_npz(path, **arrays):
    """Write an uncompressed numeric NPZ atomically, rejecting >1 MB artifacts."""
    for value in arrays.values():
        a = np.asarray(value)
        if a.dtype.hasobject:
            raise ValueError('object arrays are not portable scientific artifacts')
    if sum(np.asarray(a).nbytes for a in arrays.values()) > MAX_TRANSIENT:
        raise ValueError('array payload exceeds the atomic temporary limit')
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    _atomic_bytes(path, buffer.getvalue())


def write_array_chunks(path, array, max_payload=900_000):
    """Atomically save row chunks; return ordered paths, shapes and SHA-256.

    path is a naming prefix, e.g. truth/train_0. Each chunk is a standalone
    numeric .npy and preserves dtype. Reconstruct with concatenate(axis=0).
    Existing chunks must match exactly; corrupted or changed truth is rejected.
    """
    a = np.asarray(array)
    if a.dtype.hasobject or a.ndim < 1:
        raise ValueError('a non-scalar numeric array is required')
    if not 0 < max_payload <= 900_000:
        raise ValueError('payload must be in (0, 900000]')
    row_bytes = int(np.prod(a.shape[1:], dtype=np.int64))*a.dtype.itemsize
    if row_bytes > max_payload:
        raise ValueError('one row exceeds chunk payload budget')
    rows = max(1, max_payload // max(row_bytes, 1))
    prefix = Path(path)
    manifest = []
    starts = range(0, len(a), rows) if len(a) else [0]
    for number, start in enumerate(starts):
        chunk = np.ascontiguousarray(a[start:start+rows])
        buffer = io.BytesIO()
        np.save(buffer, chunk, allow_pickle=False)
        content = buffer.getvalue()
        digest = hashlib.sha256(content).hexdigest()
        target = prefix.parent / (prefix.name+f'.{number:04d}.npy')
        if target.exists():
            if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise RuntimeError('existing truth chunk differs; preserve and investigate')
        else:
            _atomic_bytes(target, content)
        manifest.append(dict(path=str(target), shape=list(chunk.shape),
                             dtype=chunk.dtype.str, bytes=len(content), sha256=digest))
    # A shortened or differently chunked resumed dataset must not leave stale
    # tails silently attached to its truth series.
    expected = {Path(x['path']).name for x in manifest}
    if any(p.name not in expected for p in prefix.parent.glob(prefix.name+'.*.npy')):
        raise RuntimeError('unexpected existing truth chunks; preserve and investigate')
    return manifest


def validate_condition(arrays, idx_prefix, expected_hash=None):
    """Validate and hash one completed condition without a full campaign copy.

    Conditions retain regime/object/repetition/method axes. Hash includes field
    names, shapes, dtype and C-order bytes. Missing status (128), unknown flags,
    finite failed metrics or nonfinite/negative successful metrics are rejected.
    Nonzero completed status bits are 1 numerical, 2 metric, 4 family failure.
    """
    if any(name not in arrays for name in METRICS+('status',)):
        raise ValueError('missing metric or status array')
    index = tuple(idx_prefix)
    status = np.asarray(arrays['status'][index])
    if status.dtype != np.dtype('uint8') or status.size == 0:
        raise ValueError('nonempty uint8 status required')
    if np.any(status & np.uint8(248)):
        raise ValueError('missing, corrupt or unknown condition status')
    digest = hashlib.sha256()
    for name in METRICS+('status',):
        data = np.asarray(arrays[name][index])
        if data.shape != status.shape:
            raise ValueError('condition metric/status shapes differ')
        if name != 'status' and data.dtype != np.dtype('float64'):
            raise ValueError('metrics must be float64')
        digest.update(name.encode('ascii')+b'\0')
        digest.update(repr(data.shape).encode('ascii')+b'\0'+data.dtype.str.encode('ascii')+b'\0')
        # Iterating the first remaining axis bounds a copy to one regime, even
        # if the original memmap slice is not contiguous.
        slabs = zip(data, status) if data.ndim else [(data,status)]
        for slab, flags in slabs:
            if name != 'status':
                good = flags == 0
                if np.any(~np.isfinite(slab[good])) or np.any(slab[good] < 0):
                    raise ValueError('invalid metric marked successful')
                if np.any(slab[~good] != np.inf):
                    raise ValueError('failed rows must retain positive infinite metrics')
            digest.update(np.ascontiguousarray(slab).tobytes())
    result = digest.hexdigest()
    if expected_hash is not None and result != expected_hash:
        raise ValueError('completed condition SHA-256 mismatch')
    return result

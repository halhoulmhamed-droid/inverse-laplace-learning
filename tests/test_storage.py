"""Integrity and bounded writes, including interrupted replacement."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from laplace_learning.storage import (atomic_npz, bounded_disk_usage,
    write_array_chunks, validate_condition, METRICS)


class StorageTests(unittest.TestCase):
    def test_atomic_interruption_keeps_original(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'a.npz'
            atomic_npz(path, x=np.arange(4))
            old = path.read_bytes()
            with patch('laplace_learning.storage.os.replace', side_effect=OSError('interrupted')):
                with self.assertRaises(OSError):
                    atomic_npz(path, x=np.arange(5))
            self.assertEqual(old, path.read_bytes())
            self.assertEqual(list(Path(folder).glob('*.tmp')), [])
            with self.assertRaises(ValueError):
                atomic_npz(path, x=np.zeros(130000))
            self.assertEqual(old, path.read_bytes())

    def test_chunks_and_corruption(self):
        with tempfile.TemporaryDirectory() as folder:
            values = np.arange(1000*128,dtype=float).reshape(1000,128)
            prefix = Path(folder)/'truth'
            records = write_array_chunks(prefix,values)
            self.assertEqual(len(records),2)
            self.assertTrue(all(r['bytes'] < 1_000_000 for r in records))
            restored = np.concatenate([np.load(r['path']) for r in records])
            np.testing.assert_array_equal(restored,values)
            self.assertEqual(records,write_array_chunks(prefix,values))
            Path(records[0]['path']).write_bytes(b'corruption')
            with self.assertRaises(RuntimeError):
                write_array_chunks(prefix,values)

    def test_shared_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder); (p/'pilot').mkdir(); (p/'campaign').mkdir()
            (p/'pilot'/'a').write_bytes(b'123')
            (p/'campaign'/'b').write_bytes(b'4567')
            self.assertEqual(bounded_disk_usage(p,planned_bytes=3,limit=10),7)
            with self.assertRaises(RuntimeError):
                bounded_disk_usage(p,planned_bytes=4,limit=10)

    def test_condition_hash_and_semantic_corruption(self):
        shape=(1,1,1,2,3,2,11)
        arrays={n:np.ones(shape) for n in METRICS}
        arrays['status']=np.zeros(shape,dtype=np.uint8)
        digest=validate_condition(arrays,(0,0,0))
        self.assertEqual(digest,validate_condition(arrays,(0,0,0),digest))
        arrays['visible'][0,0,0,0,0,0,0]=2
        with self.assertRaises(ValueError):
            validate_condition(arrays,(0,0,0),digest)
        arrays['visible'][0,0,0,0,0,0,0]=np.nan
        with self.assertRaises(ValueError): validate_condition(arrays,(0,0,0))
        arrays['status'][0,0,0,0,0,0,0]=1
        for n in METRICS: arrays[n][0,0,0,0,0,0,0]=np.inf
        validate_condition(arrays,(0,0,0))
        arrays['status'][0,0,0,0,0,0,0]=128
        with self.assertRaises(ValueError): validate_condition(arrays,(0,0,0))


if __name__=='__main__':
    unittest.main()

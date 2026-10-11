import unittest
import numpy as np
from laplace_learning.data import generator, profiles, noise, design, datasets


class DataTests(unittest.TestCase):
    def test_pilot_partition_ids(self):
        cfg = design("pilot")
        train, val, nominal, rough, vi, ti = datasets(cfg)
        self.assertEqual(vi, list(range(24, 32)))
        self.assertEqual(ti, list(range(32, 40)))
        self.assertEqual(len(set(range(24)) | set(vi) | set(ti)
                             | set(range(40, 56))), 56)
        np.testing.assert_array_equal(val[0], profiles(0, 0, [24], 64)[0])
        self.assertFalse(np.array_equal(train[0][0], val[0]))

    def test_seed_reproduction_and_independence(self):
        a = noise(2, 0, [0, 1], [0, 1])
        np.testing.assert_array_equal(a, noise(2, 0, [0, 1], [0, 1]))
        self.assertFalse(np.array_equal(a[:, 0], a[:, 1]))
        self.assertFalse(np.array_equal(a, noise(3, 0, [0, 1], [0, 1])))

    def test_no_profile_normalization(self):
        c = profiles(1, 0, range(10), 128)
        self.assertGreater(np.std(np.linalg.norm(c, axis=1)), 0.1)
        with self.assertRaises(ValueError):
            generator(-1, 0, 0)

    def test_counts(self):
        for mode, expected in [("full", 2376000), ("fallback", 267300)]:
            c = design(mode)
            self.assertEqual(c["games"] * len(c["masks"]) * len(c["levels"])
                             * len(c["regimes"]) * c["test"]
                             * c["repetitions"] * 11, expected)

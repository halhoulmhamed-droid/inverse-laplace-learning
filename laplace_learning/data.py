"""Disjoint physical objects and reproducible, paired sensor noise."""
import numpy as np

ROOT_SEED = 20261011


def generator(partition, game, obj, channel=0, repetition=0):
    """The six-integer seed tuple is part of the experimental protocol."""
    values = [ROOT_SEED, partition, game, obj, channel, repetition]
    if any(int(x) != x or x < 0 for x in values):
        raise ValueError("Seed components must be nonnegative integers")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(values)))


def profiles(partition, game, ids, K, p=2.5):
    weights = (1.0 + np.arange(K, dtype=float) ** 2) ** (-p / 2)
    weights /= np.sqrt(np.sum(weights ** 2))
    return np.asarray([generator(partition, game, int(i)).standard_normal(K)
                       * weights for i in ids])


def noise(partition, game, ids, repetitions, sensors=64):
    return np.asarray([[generator(partition, game, int(i), 1, int(r))
                        .standard_normal(sensors) for r in repetitions]
                       for i in ids])


def design(mode):
    if mode == "full":
        return dict(mode=mode, n=64, K=128, nx=64, ny=16, train=1000,
                    validation=200, test=200, repetitions=10, games=3,
                    masks=["full", "even", "contiguous"],
                    levels=[0.001, 0.01, 0.05],
                    regimes=["nominal", "rough", "height_2", "height_5"])
    if mode == "fallback":
        return dict(mode=mode, n=32, K=64, nx=32, ny=8, train=300,
                    validation=60, test=60, repetitions=5, games=3,
                    masks=["full", "even", "contiguous"],
                    levels=[0.001, 0.01, 0.05],
                    regimes=["nominal", "rough", "height_5"])
    if mode == "pilot":
        return dict(mode=mode, n=32, K=64, nx=32, ny=8, train=24,
                    validation=8, test=8, repetitions=2, games=1,
                    masks=["full", "contiguous"], levels=[0.01],
                    regimes=["nominal", "height_5"])
    raise ValueError("Unknown design")


def datasets(config):
    K = config["K"]
    pilot = config["mode"] == "pilot"
    train = [profiles(0 if pilot else 1, g,
                      range(config["train"]), K)
             for g in range(config["games"])]
    vi = range(24, 32) if pilot else range(config["validation"])
    ti = range(32, 40) if pilot else range(config["test"])
    val = profiles(0 if pilot else 2, 0, vi, K)
    nominal = profiles(0 if pilot else 3, 0, ti, K)
    rough = None if pilot else profiles(4, 0, ti, K, p=1.5)
    return train, val, nominal, rough, list(vi), list(ti)

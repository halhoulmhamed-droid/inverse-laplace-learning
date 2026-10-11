# Laplace data completion by learned spectral filtering

**Projet scientifique de Master — Halhoul Mhamed**<br>
Université Abdelmalek Essaâdi · Faculté Polydisciplinaire de Larache<br>
Master « Modélisation mathématique et applications et apprentissages » · M2 en cours<br>
**Professeur encadrant : Pr. Chakir Tajani**

This project compares eleven affine linear estimators for reconstructing the upper boundary trace of a harmonic field from noisy lower-boundary averages. The lower and lateral fluxes are assumed to be zero. Training estimates spectral gains and covariances; validation uses two independent observations, never validation truth. No neural network is required.

The analytic data generator and finite-difference inverse are independent. Experiments include partial sensor coverage, rougher test traces, and height mismatch. The complete campaign contains **2,376,000 reconstructions** across three independent training sets, with no numerical failures. The predeclared 10% reduction in the 95th error percentile from capping learned signed gains **was not established**. These synthetic results are not an algorithmic novelty claim or a physical validation.

## Scientific material

- [Mathematical model and proofs](docs/theory.md)
- [Pre-experiment protocol](docs/protocol.md)
- [Sources and reading limitations](docs/sources.md)
- [Author's preparation note, in French](docs/defense.md)
- [Primary comparisons and bootstrap conclusions](results/analysis/primary.json)
- [All comparison summaries](results/analysis/summary.csv)
- [Measured campaign checks and costs](results/full/report.json)

The French report, slides and figure provenance are prepared from these archived results. Their final release status is recorded with the document delivery.

## Reproduce

Python 3.12, NumPy, SciPy and Matplotlib suffice. Install the pinned free dependencies in an isolated environment. No FreeFem, GPU or paid service is needed.

```bash
python -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/check_tests.py
python -m laplace_learning.experiment --mode pilot --output results/reproduction_pilot --tests-report results/tests.json
```

For a fresh full reproduction, use a **new empty output directory** and a successful pilot report:

```bash
python -m laplace_learning.experiment --mode full --output results/reproduction_full --tests-report results/tests.json --pilot-report results/reproduction_pilot/report.json
python scripts/analyze.py --help
python scripts/figures.py --help
```

The 200 MB retained-output limit covers the common results parent directory. Keeping the archived campaign and an additional full copy there exceeds the preflight allowance. For an independent full reproduction, use an empty separate experiment directory, retaining the published archive unchanged. Configuration or source-hash mismatches intentionally reject resuming an old run. Completed conditions are checked by SHA-256 and semantic status checks; unfinished conditions can resume within the remaining budget.

The seeds use NumPy PCG64 and explicit object partitions. Truth chunks concatenate in manifest order along axis 0. Metric arrays use the axes recorded in `configuration.json`. Status 0 is success; bits 1, 2 and 4 mark numerical, metric and family failures; 128 means missing. Nonfinite failed reconstructions remain in the statistics as infinite errors.

Measured timings reflect this small linear algebra workload and this environment, not a general runtime promise. The recorded cumulative budget includes conservative provisions. Bootstrap uncertainty is conditional on each training set; it is not a simultaneous confidence statement over all three training sets.

## Attribution and use

Third-party teaching materials and thesis PDFs are not redistributed. References retain their authorship and reading limitations. No open-source licence has yet been selected for this repository; do not infer a licence from public availability. Institutional marks, when provided with verified provenance, remain the property of their institutions. No institutional endorsement or supervisor approval is claimed.

## Contact

halhoul.mhamed@etu.uae.ac.ma<br>
halhoulmhamed@gmail.com<br>
https://github.com/halhoulmhamed-droid/inverse-laplace-learning

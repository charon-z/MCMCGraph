# Computational-performance benchmark

This directory contains the reproducible timing and peak-memory benchmark used
to characterize BPFC software requirements. It compares BPFC-MCMC with the
same-model BPFC-EM optimizer, a full-covariance Gaussian mixture (`mclust`),
K-means, and Ward agglomerative clustering.

All methods are timed in fresh R processes on the same machine. Classical
comparators include Legendre-feature extraction in their end-to-end timing.
The shell runner fixes common numerical-library thread variables to one. GNU
`/usr/bin/time -v` records peak resident memory independently for every run.

From an isolated benchmark directory containing `source/BPFC`, `rlib`, `logs`,
`results`, and `data`:

```bash
export BPFC_BENCH_ROOT="$PWD"
export BPFC_BENCH_SOURCE="$PWD/source/BPFC"
export BPFC_BENCH_RLIB="$PWD/rlib"
export BPFC_BENCH_MOUSE_INPUT="$PWD/data/individual_combined.txt"
bash source/BPFC/inst/benchmark/run_benchmark_suite.sh \
  > logs/suite.log 2>&1
```

The runner is resumable: a successful per-run CSV and status file cause that
configuration to be skipped on restart. The mouse analysis is attempted only
when `BPFC_BENCH_MOUSE_INPUT` exists. Raw and summarized tables are written to
`results/computational_performance_raw.csv` and
`results/computational_performance_summary.csv`.

The 30,000-iteration mouse fit is an end-to-end production benchmark. Do not
replace it with a linearly extrapolated value unless the manuscript labels the
result explicitly as an estimate.

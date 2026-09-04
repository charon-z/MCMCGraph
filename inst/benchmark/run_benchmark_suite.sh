#!/usr/bin/env bash
set -u

BENCH_ROOT="${BPFC_BENCH_ROOT:-$(pwd)}"
SOURCE_DIR="${BPFC_BENCH_SOURCE:-${BENCH_ROOT}/source/BPFC}"
R_LIBRARY="${BPFC_BENCH_RLIB:-${BENCH_ROOT}/rlib}"
MOUSE_INPUT="${BPFC_BENCH_MOUSE_INPUT:-${BENCH_ROOT}/data/individual_combined.txt}"
RUNNER="${SOURCE_DIR}/inst/benchmark/benchmark_one.R"

mkdir -p "${BENCH_ROOT}/logs" "${BENCH_ROOT}/results/runs"

export R_LIBS_USER="${R_LIBRARY}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

{
  date --iso-8601=seconds
  uname -a
  lscpu
  free -h
  Rscript --version
  git -C "${SOURCE_DIR}" rev-parse HEAD
} > "${BENCH_ROOT}/results/system_info.txt" 2>&1

run_one() {
  method="$1"
  dataset="$2"
  n="$3"
  J="$4"
  niter="$5"
  repeat_id="$6"
  seed="$7"
  run_id="${dataset}_n${n}_J${J}_${method}_iter${niter}_rep${repeat_id}"
  csv_file="${BENCH_ROOT}/results/runs/${run_id}.csv"
  time_file="${BENCH_ROOT}/results/runs/${run_id}.time.txt"
  log_file="${BENCH_ROOT}/logs/${run_id}.log"
  status_file="${BENCH_ROOT}/results/runs/${run_id}.status"

  if [[ -s "${csv_file}" && -s "${status_file}" ]] &&
     [[ "$(cat "${status_file}")" == "0" ]]; then
    echo "SKIP completed ${run_id}"
    return 0
  fi

  echo "START ${run_id} $(date --iso-8601=seconds)"
  /usr/bin/time -v -o "${time_file}" \
    env \
      BPFC_BENCH_METHOD="${method}" \
      BPFC_BENCH_DATASET="${dataset}" \
      BPFC_BENCH_N="${n}" \
      BPFC_BENCH_J="${J}" \
      BPFC_BENCH_NITER="${niter}" \
      BPFC_BENCH_REPEAT="${repeat_id}" \
      BPFC_BENCH_SEED="${seed}" \
      BPFC_BENCH_OUTPUT="${csv_file}" \
      BPFC_BENCH_MOUSE_INPUT="${MOUSE_INPUT}" \
      BPFC_BENCH_SOURCE="${SOURCE_DIR}" \
      Rscript "${RUNNER}" > "${log_file}" 2>&1
  status=$?
  printf '%s\n' "${status}" > "${status_file}"
  echo "END ${run_id} status=${status} $(date --iso-8601=seconds)"
}

# Method comparison at the principal simulated setting. Three repeats provide
# stable median timings without repeating the full scientific simulation.
for repeat_id in 1 2 3; do
  seed=$((20260904 + repeat_id))
  for method in bpfc bpfc_em gmm kmeans agglomerative; do
    run_one "${method}" sim 3000 8 5000 "${repeat_id}" "${seed}"
  done
done

# BPFC scaling with feature number and component number. The n=3000, J=8 cell
# above is shared with the method comparison.
run_one bpfc sim 1000 8 5000 1 20261001
run_one bpfc sim 10000 8 5000 1 20261002
run_one bpfc sim 3000 5 5000 1 20261003
run_one bpfc sim 3000 15 5000 1 20261004

# Real-data production configuration. Run once because this 30,000-iteration
# fit is intentionally the expensive reproduction benchmark.
if [[ -f "${MOUSE_INPUT}" ]]; then
  run_one bpfc mouse 11833 15 30000 1 20261101
  for method in bpfc_em gmm kmeans agglomerative; do
    run_one "${method}" mouse 11833 15 0 1 20261101
  done
else
  echo "WAITING: mouse input not found at ${MOUSE_INPUT}"
fi

Rscript "${SOURCE_DIR}/inst/benchmark/summarize_benchmark.R" \
  "${BENCH_ROOT}/results/runs" \
  "${BENCH_ROOT}/results"

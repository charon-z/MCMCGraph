#!/usr/bin/env bash
# run_full.sh -------------------------------------------------------------
# Launches the full Task-1 simulation study as a queue of checkpointed jobs
# with bounded parallelism. Each (experiment, scenario, replicate) cell writes
# its own .rds and is skipped if present, so this script is safe to re-run and
# resume after interruption.
#
# Usage:   bash run_full.sh [N_PARALLEL]      (default 6)
# Env:     MCG_ROOT (default: parent of this script's ../../..)
#          MCG_R    (default 50)              total replicates
#
# Scale is tuned to finish overnight on a 16 GB machine without swapping:
# headline experiments (1A, 1B, 1D) use niter = 5000-6000; the supplementary 1C
# ablation uses niter = 3000 and R = 8. For the as-published R = 50 / niter =
# 8000 setting, raise these literals (a larger-RAM host or server is advised).

set -u
NP="${1:-6}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MCG_ROOT="${MCG_ROOT:-$(cd "$HERE/../../.." && pwd)}"
export MCG_LIB_DIR="$HERE"
R_TOTAL="${MCG_R:-50}"
LOGDIR="$MCG_ROOT/results/reproduce/logs"; mkdir -p "$LOGDIR"
RS="$HERE/01_run_simulations.R"

half1="1:$((R_TOTAL/2))"; half2="$((R_TOTAL/2+1)):$R_TOTAL"

# Build the job list: each line = "TAG|EXP|SCENARIOS|REPS|NITER|R"
JOBS=()
for K in 3 5 8; do
  JOBS+=("1A_K${K}_a|1A|c($K)|$half1|5000|$R_TOTAL")
  JOBS+=("1A_K${K}_b|1A|c($K)|$half2|5000|$R_TOTAL")
  JOBS+=("1Ainit_K${K}|1Ainit|c($K)|1|5000|$R_TOTAL")
  JOBS+=("1Bmain_K${K}_a|1Bmain|c($K)|$half1|5000|$R_TOTAL")
  JOBS+=("1Bmain_K${K}_b|1Bmain|c($K)|$half2|5000|$R_TOTAL")
done
JOBS+=("1Bmirror_a|1Bmirror|c(0)|$half1|5000|$R_TOTAL")
JOBS+=("1Bmirror_b|1Bmirror|c(0)|$half2|5000|$R_TOTAL")
for K in 3 5 8; do
  JOBS+=("1C_K${K}|1C|c($K)|1:8|3000|8")
done
JOBS+=("1D|1D|c(3,5,8)|1|6000|3")

echo "[run_full] $(date)  ROOT=$MCG_ROOT  parallel=$NP  jobs=${#JOBS[@]}"

run_one() {
  local spec="$1"
  IFS='|' read -r tag exp scn reps niter rr <<< "$spec"
  local out="$LOGDIR/${tag}.log"
  if [ -f "$LOGDIR/${tag}.done" ]; then echo "[skip] $tag"; return 0; fi
  echo "[start] $tag  exp=$exp scn=$scn reps=$reps niter=$niter"
  MCG_EXP="$exp" MCG_SCENARIOS="$scn" MCG_REPS="$reps" MCG_NITER="$niter" MCG_R="$rr" \
    Rscript "$RS" > "$out" 2>&1 && touch "$LOGDIR/${tag}.done"
  echo "[end]   $tag exit=$?"
}
export -f run_one
export LOGDIR RS

printf "%s\n" "${JOBS[@]}" | xargs -I{} -P "$NP" bash -c 'run_one "$@"' _ {}

echo "[run_full] all jobs dispatched/complete $(date)"

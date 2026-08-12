#!/usr/bin/env bash
# server_bootstrap.sh -----------------------------------------------------
# Disconnect-proof launcher for the full BPFC Task-1 study on a remote
# host. Run this ON THE SERVER (e.g. rqzhao@222.28.118.18:~/workdir/MCMC).
# It detaches with setsid+nohup so it keeps running after you log out, and is
# checkpointed/resumable, so re-running just continues where it stopped.
#
# Usage (on the server, inside the unpacked bundle directory):
#   bash server_bootstrap.sh            # install deps + launch full run detached
#   bash server_bootstrap.sh status     # show progress
#   bash server_bootstrap.sh tables     # build tables+figures from whatever is done
#
# Scale: full as-published R = 50 / niter = 8000 (override with MCG_R / NP env).

set -u
ACTION="${1:-run}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# bundle layout: <root>/BPFC/...  and  <root>/benchmark/data/sim_truth_*.rds
ROOT="$(cd "$HERE/../../.." && pwd)"
export MCG_ROOT="$ROOT"
export MCG_LIB_DIR="$HERE"
NP="${NP:-6}"                      # parallel workers (servers usually have RAM)
export MCG_R="${MCG_R:-50}"
LOGDIR="$ROOT/results/reproduce/logs"; mkdir -p "$LOGDIR"
MASTER="$ROOT/results/reproduce/server_master.log"

deps() {
  Rscript -e '
    options(repos=c(CRAN="https://cloud.r-project.org"))
    need <- c("nimble","mclust","clue","e1071","xgboost","coda","ggplot2","igraph")
    for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
    cat("deps:\n"); print(sapply(need, function(p) requireNamespace(p, quietly=TRUE)))
  '
  R CMD INSTALL "$ROOT/BPFC"
}

case "$ACTION" in
  status)
    echo "cells: 1A=$(ls "$ROOT"/results/reproduce/sim1A/K*_rep*.rds 2>/dev/null|wc -l) initsens=$(ls "$ROOT"/results/reproduce/sim1A/initsens*.rds 2>/dev/null|wc -l) 1B=$(ls "$ROOT"/results/reproduce/sim1B/K*_rep*.rds 2>/dev/null|wc -l) mirror=$(ls "$ROOT"/results/reproduce/sim1B/mirror*.rds 2>/dev/null|wc -l) 1C=$(ls "$ROOT"/results/reproduce/sim1C/*.rds 2>/dev/null|wc -l) 1D=$(ls "$ROOT"/results/reproduce/sim1D/*.rds 2>/dev/null|wc -l)"
    echo "done_jobs=$(ls "$LOGDIR"/*.done 2>/dev/null|wc -l)/18  R_procs=$(pgrep -f 01_run_simulations|wc -l)"
    tail -2 "$MASTER" 2>/dev/null
    ;;
  tables)
    Rscript "$HERE/02_make_tables.R"; Rscript "$HERE/03_make_figures.R"
    ;;
  deps) deps ;;
  run)
    echo "[bootstrap] installing deps ($(date))"; deps
    echo "[bootstrap] launching detached full run, NP=$NP R=$MCG_R ($(date))"
    # detached, survives logout; resumable via .done markers
    MCG_R="$MCG_R" setsid nohup bash "$HERE/run_full.sh" "$NP" > "$MASTER" 2>&1 < /dev/null &
    echo "[bootstrap] PID $! ; tail -f $MASTER"
    ;;
  *) echo "unknown action: $ACTION"; exit 1 ;;
esac

#!/usr/bin/env Rscript
# 02_make_tables.R --------------------------------------------------------
# Aggregates the per-cell result files written by 01_run_simulations.R into the
# paper tables: mean +/- SD (and median [IQR]) per scenario x method x metric,
# the joint/concat/independent ablation table, the secondary-ablation tables,
# and the MCMC-vs-EM paired Wilcoxon signed-rank tests.
#
# Every table is wrapped so a missing/partial input never aborts the others.
#
# Outputs (results/reproduce/tables/):
#   table_1A_methods.csv          mean+/-SD + median[IQR], 5 methods x 3 scenarios
#   table_1A_wilcoxon.csv         MCMC vs EM paired Wilcoxon p per metric/scenario
#   table_1A_initsens.csv         EM vs MCMC ARI spread across inits
#   table_1B_ablation.csv         joint/concat/independent (+ mirror) mean+/-SD
#   table_1B_wilcoxon.csv         joint vs concat / joint vs independent paired p
#   table_1C_covariance.csv, table_1C_order.csv, table_1C_prior.csv
#   table_1D_convergence.csv      Rhat / ESS summaries

src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "MCMCGraph/inst/reproduce"
source(file.path(src_dir, "config.R"))

metric_names <- c("ACC", "ARI", "NMI", "Macro_F1", "Balanced_ACC", "SmallCluster_Recall")
TAB <- file.path(MCG_OUT_DIR, "tables")
try_write <- function(expr) tryCatch(expr, error = function(e) message("  [skip] ", conditionMessage(e)))

as_scalar_numeric <- function(x) {
  if (!is.list(x)) return(as.numeric(x))
  vapply(x, function(v) {
    if (!length(v)) return(NA_real_)
    suppressWarnings(as.numeric(unlist(v, recursive = TRUE, use.names = FALSE)[1]))
  }, numeric(1))
}

normalize_cell_types <- function(x) {
  if (is.null(x) || !is.data.frame(x)) return(x)
  numeric_cols <- intersect(
    c("scenario", "rep", "init", "seed", "K", "niter", "n_chains", metric_names),
    names(x)
  )
  for (nm in numeric_cols) x[[nm]] <- as_scalar_numeric(x[[nm]])
  x
}

read_cells <- function(sub, pattern) {
  files <- list.files(file.path(MCG_OUT_DIR, sub), pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  normalize_cell_types(do.call(rbind, lapply(files, function(f)
    tryCatch(readRDS(f), error = function(e) NULL))))
}

scen_label <- function(s) {
  lab <- c(`3` = "K3 (balanced)", `5` = "K5 (unbalanced/overlap)",
           `8` = "K8 (small clusters)", `0` = "Mirror (orthogonal)",
           `90` = "Mirror (pairing)")[as.character(s)]
  ifelse(is.na(lab), paste0("scenario ", s), lab)
}
fmt_ms   <- function(x) sprintf("%.3f +/- %.3f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_miqr <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x, na.rm = TRUE),
                                quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))

summ_table <- function(df, group_col = "method") {
  out <- list()
  for (s in sort(unique(df$scenario))) {
    sub <- df[df$scenario == s, , drop = FALSE]
    for (g in unique(sub[[group_col]])) {
      gg <- sub[sub[[group_col]] == g, , drop = FALSE]
      r <- data.frame(scenario = scen_label(s), group = g, n_rep = nrow(gg),
                      stringsAsFactors = FALSE)
      for (m in metric_names) {
        r[[paste0(m, "_meanSD")]] <- fmt_ms(gg[[m]])
        r[[paste0(m, "_medIQR")]] <- fmt_miqr(gg[[m]])
      }
      out[[length(out) + 1]] <- r
    }
  }
  do.call(rbind, out)
}

paired_wilcox <- function(df, a, b, group_col = "method", by = "rep") {
  out <- list()
  for (s in sort(unique(df$scenario))) {
    sub <- df[df$scenario == s, , drop = FALSE]
    A <- sub[sub[[group_col]] == a, , drop = FALSE]; B <- sub[sub[[group_col]] == b, , drop = FALSE]
    if (!nrow(A) || !nrow(B)) next
    mrg <- merge(A, B, by = by, suffixes = c("_a", "_b"))
    for (m in metric_names) {
      xa <- mrg[[paste0(m, "_a")]]; xb <- mrg[[paste0(m, "_b")]]
      ok <- is.finite(xa) & is.finite(xb)
      p <- if (sum(ok) >= 2 && any(xa[ok] != xb[ok]))
        suppressWarnings(wilcox.test(xa[ok], xb[ok], paired = TRUE)$p.value) else NA_real_
      out[[length(out) + 1]] <- data.frame(
        scenario = scen_label(s), comparison = sprintf("%s vs %s", a, b),
        metric = m, n_pair = sum(ok),
        median_diff = if (any(ok)) median(xa[ok] - xb[ok]) else NA_real_,
        p_value = p, stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

## ---------------- 1A ----------------
a1 <- read_cells("sim1A", "^K\\d+_rep")
if (!is.null(a1)) {
  main <- a1[a1$method %in% c("MCMC", "EM", "LOP-kmeans", "SVM-RBF", "XGBoost"), , drop = FALSE]
  try_write(write.csv(summ_table(main), file.path(TAB, "table_1A_methods.csv"), row.names = FALSE))
  try_write(write.csv(paired_wilcox(main, "MCMC", "EM"), file.path(TAB, "table_1A_wilcoxon.csv"), row.names = FALSE))
  ei <- a1[a1$method == "EM_init", , drop = FALSE]
  if (nrow(ei)) try_write({
    sp <- aggregate(ARI ~ scenario + rep, ei, function(x) c(sd = sd(x), rng = diff(range(x))))
    spdf <- data.frame(scenario = scen_label(sp$scenario), rep = sp$rep,
                       ARI_sd = sp$ARI[, "sd"], ARI_range = sp$ARI[, "rng"])
    write.csv(spdf, file.path(TAB, "table_1A_EM_init_spread.csv"), row.names = FALSE)
  })
  message("[1A] tables written (", nrow(main), " method rows)")
}
is1 <- read_cells("sim1A", "^initsens")
if (!is.null(is1)) try_write({
  agg <- aggregate(ARI ~ scenario + method, is1, function(x)
    c(mean = mean(x), sd = sd(x), min = min(x), max = max(x)))
  out <- data.frame(scenario = scen_label(agg$scenario), method = agg$method,
                    ARI_mean = agg$ARI[, "mean"], ARI_sd = agg$ARI[, "sd"],
                    ARI_min = agg$ARI[, "min"], ARI_max = agg$ARI[, "max"])
  write.csv(out, file.path(TAB, "table_1A_initsens.csv"), row.names = FALSE)
})

## ---------------- 1B ----------------
b_all <- rbind(read_cells("sim1B", "^K\\d+_rep"), read_cells("sim1B", "^mirror"))
if (!is.null(b_all)) {
  ab <- b_all[b_all$arm %in% c("joint", "concat", "indep", "indep_crossprod"), , drop = FALSE]
  try_write(write.csv(summ_table(ab, "arm"), file.path(TAB, "table_1B_ablation.csv"), row.names = FALSE))
  try_write({
    wj <- rbind(paired_wilcox(ab, "joint", "concat", "arm"),
                paired_wilcox(ab, "joint", "indep", "arm"))
    write.csv(wj, file.path(TAB, "table_1B_wilcoxon.csv"), row.names = FALSE)
  })
  message("[1B] ablation tables written")
}

## ---------------- 1C ----------------
c1 <- read_cells("sim1C", "^K\\d+_rep")
if (!is.null(c1)) {
  for (arm in c("cov", "order", "prior")) {
    sub <- c1[c1$arm == arm, , drop = FALSE]
    if (nrow(sub)) try_write(write.csv(summ_table(sub, "method"),
                             file.path(TAB, sprintf("table_1C_%s.csv", arm)), row.names = FALSE))
  }
  message("[1C] secondary-ablation tables written")
}

## ---------------- 1D ----------------
d_files <- list.files(file.path(MCG_OUT_DIR, "sim1D"), pattern = "^conv", full.names = TRUE)
if (length(d_files)) try_write({
  rh_of <- function(g) if (!is.null(g)) g$psrf[, 1] else NA_real_
  rows <- lapply(d_files, function(f) {
    x   <- readRDS(f)
    inv <- rh_of(x$gelman_inv)             # label-invariant params (phi, v_sq, v_raw)
    aln <- rh_of(x$gelman_aligned)         # component params after label alignment
    raw <- rh_of(x$gelman)                 # all params, no alignment (transparency)
    phi <- x$phi_by_chain
    data.frame(
      K = x$K, n_chains = x$n_chains, niter = x$niter,
      # primary convergence target: identifiable (label-invariant) parameters
      n_inv = length(inv),
      Rhat_inv_max = max(inv, na.rm = TRUE), Rhat_inv_median = median(inv, na.rm = TRUE),
      Rhat_inv_frac_below_1.01 = mean(inv < 1.01, na.rm = TRUE),
      ESS_inv_min = if (!is.null(x$ess_inv)) min(x$ess_inv, na.rm = TRUE) else NA,
      # component params after between-chain label alignment
      Rhat_aligned_max = max(aln, na.rm = TRUE), Rhat_aligned_median = median(aln, na.rm = TRUE),
      ESS_aligned_min = if (!is.null(x$ess_aligned)) min(x$ess_aligned, na.rm = TRUE) else NA,
      # multimodality flag: spread of per-chain posterior-mean phi (shared scalar)
      phi_chain_spread = if (length(phi) > 1) max(phi) - min(phi) else NA,
      # raw (no alignment) kept for transparency
      Rhat_raw_max = max(raw, na.rm = TRUE))
  })
  write.csv(do.call(rbind, rows), file.path(TAB, "table_1D_convergence.csv"), row.names = FALSE)
  message("[1D] convergence table written")
})
message("[02] tables in ", TAB)

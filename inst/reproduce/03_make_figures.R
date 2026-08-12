#!/usr/bin/env Rscript
# 03_make_figures.R -------------------------------------------------------
# Builds the figures for Task 1 from the per-cell result files and (for the
# mechanism panels) freshly generated example data.
#
# Outputs (results/reproduce/figures/):
#   fig_1A_boxplot_ARI.{pdf,png}        ARI by method, faceted by scenario
#   fig_1A_init_sensitivity.{pdf,png}   EM vs MCMC ARI scatter across inits
#   fig_1B_ablation_ARI.{pdf,png}       joint/concat/independent boxplots
#   fig_1B_mirror_mechanism.{pdf,png}   trajectory panels explaining the mirror
#   fig_1D_trace.{pdf,png}              representative trace plots
#
# Uses base graphics + ggplot2 if available; degrades gracefully to base.

src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "BPFC/inst/reproduce"
source(file.path(src_dir, "config.R"))
source(file.path(src_dir, "lib_sim.R"))
has_gg <- requireNamespace("ggplot2", quietly = TRUE)
FIG <- file.path(MCG_OUT_DIR, "figures")
metric_names <- c("ACC", "ARI", "NMI", "Macro_F1", "Balanced_ACC", "SmallCluster_Recall")

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
  normalize_cell_types(do.call(rbind, lapply(files, readRDS)))
}
save_both <- function(name, expr, w = 8, h = 5) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); print(expr); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w * 120, height = h * 120, res = 120); print(expr); dev.off()
}
scen_lab <- function(s) {
  m <- c(`3` = "K3 balanced", `5` = "K5 unbalanced", `8` = "K8 small clusters",
         `0` = "Mirror (orth.)", `90` = "Mirror (pairing)")
  factor(m[as.character(s)], levels = unname(m))
}

## ---- Fig 1A: ARI boxplot by method x scenario ----
a1 <- read_cells("sim1A", "^K\\d+_rep")
if (!is.null(a1) && has_gg) {
  library(ggplot2)
  d <- a1[a1$method %in% c("MCMC","EM","LOP-kmeans","SVM-RBF","XGBoost"), ]
  d$scenario <- scen_lab(d$scenario)
  d$method <- factor(d$method, levels = c("MCMC","EM","LOP-kmeans","SVM-RBF","XGBoost"))
  p <- ggplot(d, aes(method, ARI, fill = method)) +
    geom_boxplot(outlier.size = 0.6, width = 0.65) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Clustering accuracy across replicates") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  save_both("fig_1A_boxplot_ARI", p, w = 10, h = 4.2)
  message("[fig] 1A boxplot done")
}

## ---- Fig 1A: EM vs MCMC init sensitivity ----
is1 <- read_cells("sim1A", "^initsens")
if (!is.null(is1) && has_gg) {
  library(ggplot2)
  is1$scenario <- scen_lab(is1$scenario)
  p <- ggplot(is1, aes(method, ARI, colour = method)) +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "black", linewidth = 0.3) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Sensitivity to initialisation (one dataset, multiple starts)") +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  save_both("fig_1A_init_sensitivity", p, w = 9, h = 4)
  message("[fig] 1A init-sensitivity done")
}

## ---- Fig 1B: ablation boxplot ----
b_all <- rbind(read_cells("sim1B", "^K\\d+_rep"), read_cells("sim1B", "^mirror"))
if (!is.null(b_all) && has_gg) {
  library(ggplot2)
  d <- b_all[b_all$arm %in% c("joint","concat","indep"), ]
  d$scenario <- scen_lab(d$scenario)
  d$arm <- factor(d$arm, levels = c("joint","concat","indep"),
                  labels = c("Joint","Concatenation","Independent"))
  p <- ggplot(d, aes(arm, ARI, fill = arm)) +
    geom_boxplot(width = 0.6, outlier.size = 0.6) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Clustering-unit ablation (same Bayesian SAD-MCMC kernel)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  save_both("fig_1B_ablation_ARI", p, w = 10, h = 4.2)
  message("[fig] 1B ablation done")
}

## ---- Fig 1B: mirror mechanism (trajectory means per cluster) ----
mir <- mcg_gen_mirror(seed = mcg_seed(99L, 1L), n = 600L)
{
  mu <- mir$mu
  cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    for (ind in 1:2) {
      idx <- ((ind - 1) * 16 + 1):(ind * 16)
      matplot(1:16, t(mu[, idx]), type = "l", lty = 1, lwd = 2.5, col = cols,
              xlab = "time", ylab = "effect",
              main = sprintf("Individual %d (marginal mean curves)", ind),
              ylim = range(mu))
      if (ind == 1) legend("topleft", legend = paste("cluster", 1:4),
                           col = cols, lwd = 2.5, bty = "n", cex = 0.9)
    }
  }
  pdf(file.path(FIG, "fig_1B_mirror_mechanism.pdf"), width = 9, height = 4); draw(); dev.off()
  png(file.path(FIG, "fig_1B_mirror_mechanism.png"), width = 1080, height = 480, res = 120); draw(); dev.off()
  message("[fig] 1B mirror mechanism done")
}

## ---- Fig 1D: trace plots ----
d_files <- list.files(file.path(MCG_OUT_DIR, "sim1D"), pattern = "^conv", full.names = TRUE)
if (length(d_files)) {
  x <- readRDS(d_files[1])
  ml <- x$chains
  pars <- intersect(c("phi", "phi_free", "p[1]", "p[2]", "beta[1, 1]"), coda::varnames(ml))
  if (length(pars)) {
    draw <- function() {
      par(mfrow = c(length(pars), 1), mar = c(3, 4, 2, 1))
      for (pn in pars) {
        rng <- range(sapply(ml, function(c) c[, pn]))
        plot(NA, xlim = c(1, nrow(ml[[1]])), ylim = rng, xlab = "iter", ylab = pn,
             main = paste("trace:", pn))
        for (ci in seq_along(ml)) lines(ml[[ci]][, pn], col = ci)
      }
    }
    pdf(file.path(FIG, "fig_1D_trace.pdf"), width = 7, height = 8); draw(); dev.off()
    png(file.path(FIG, "fig_1D_trace.png"), width = 840, height = 960, res = 120); draw(); dev.off()
    message("[fig] 1D trace done")
  }
}
message("[03] figures in ", FIG)

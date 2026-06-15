#!/usr/bin/env Rscript
# 03_make_figures_local.R --------------------------------------------------
# LOCAL (macOS) re-render of the Task-1 figures with Arial, large fonts
# (base_size 26) and 300-dpi sharp PNGs, so panels stay legible after
# assembly into composite figures. Mirrors 03_make_figures.R logic.

FONT <- "Arial"; BASE <- 26
src_dir <- Sys.getenv("MCG_LIB_DIR")
source(file.path(src_dir, "config.R"))
source(file.path(src_dir, "lib_sim.R"))
has_gg <- requireNamespace("ggplot2", quietly = TRUE)
FIG <- file.path(MCG_OUT_DIR, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
metric_names <- c("ACC","ARI","NMI","Macro_F1","Balanced_ACC","SmallCluster_Recall")

as_scalar_numeric <- function(x) {
  if (!is.list(x)) return(as.numeric(x))
  vapply(x, function(v) { if (!length(v)) return(NA_real_)
    suppressWarnings(as.numeric(unlist(v, recursive = TRUE, use.names = FALSE)[1])) }, numeric(1))
}
normalize_cell_types <- function(x) {
  if (is.null(x) || !is.data.frame(x)) return(x)
  nc <- intersect(c("scenario","rep","init","seed","K","niter","n_chains", metric_names), names(x))
  for (nm in nc) x[[nm]] <- as_scalar_numeric(x[[nm]]); x
}
read_cells <- function(sub, pattern) {
  files <- list.files(file.path(MCG_OUT_DIR, sub), pattern = pattern, full.names = TRUE)
  if (!length(files)) return(NULL)
  normalize_cell_types(do.call(rbind, lapply(files, readRDS)))
}

# sharp + REAL Arial: macOS quartz (CoreText) for BOTH pdf and 300-dpi png.
# (cairo_pdf silently substitutes Arial -> Hiragino on this fontconfig, so avoid it.)
save_both <- function(name, expr, w = 8, h = 5) {
  quartz(file = file.path(FIG, paste0(name, ".pdf")), type = "pdf", width = w, height = h, family = FONT)
  print(expr); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in",
      res = 300, type = "quartz", family = FONT)
  print(expr); dev.off()
}
save_base <- function(name, draw, w = 9, h = 4) {
  quartz(file = file.path(FIG, paste0(name, ".pdf")), type = "pdf", width = w, height = h, family = FONT); draw(); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in",
      res = 300, type = "quartz", family = FONT); draw(); dev.off()
}
scen_lab <- function(s) {
  m <- c(`3`="K3 balanced",`5`="K5 unbalanced",`8`="K8 small clusters",
         `0`="Mirror (orth.)",`90`="Mirror (pairing)")
  factor(m[as.character(s)], levels = unname(m))
}
THEME <- function() ggplot2::theme_bw(base_size = BASE, base_family = FONT)

## ---- Fig 1A: ARI boxplot by method x scenario ----
a1 <- read_cells("sim1A", "^K\\d+_rep")
if (!is.null(a1) && has_gg) {
  library(ggplot2)
  d <- a1[a1$method %in% c("MCMC","EM","LOP-kmeans","SVM-RBF","XGBoost"), ]
  d$scenario <- scen_lab(d$scenario)
  d$method <- factor(d$method, levels = c("MCMC","EM","LOP-kmeans","SVM-RBF","XGBoost"))
  p <- ggplot(d, aes(method, ARI, fill = method)) +
    geom_boxplot(outlier.size = 0.8, width = 0.65, linewidth = 0.6) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Clustering accuracy across replicates") +
    THEME() + theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  save_both("fig_1A_boxplot_ARI", p, w = 16, h = 6.5); message("[fig] 1A boxplot done")
}

## ---- Fig 1A: EM vs MCMC init sensitivity ----
is1 <- read_cells("sim1A", "^initsens")
if (!is.null(is1) && has_gg) {
  library(ggplot2)
  is1$scenario <- scen_lab(is1$scenario)
  p <- ggplot(is1, aes(method, ARI, colour = method)) +
    geom_jitter(width = 0.12, height = 0, size = 3.2, alpha = 0.85) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, colour = "black", linewidth = 0.5) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Sensitivity to initialisation") +
    THEME() + theme(legend.position = "none")
  save_both("fig_1A_init_sensitivity", p, w = 15, h = 6); message("[fig] 1A init-sensitivity done")
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
    geom_boxplot(width = 0.6, outlier.size = 0.8, linewidth = 0.6) +
    facet_wrap(~scenario, nrow = 1) +
    labs(x = NULL, y = "ARI", title = "Clustering-unit ablation") +
    THEME() + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none")
  save_both("fig_1B_ablation_ARI", p, w = 18, h = 6.5); message("[fig] 1B ablation done")
}

## ---- Fig 1B: mirror mechanism (base graphics) ----
mir <- mcg_gen_mirror(seed = mcg_seed(99L, 1L), n = 600L)
{
  mu <- mir$mu; cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a")
  draw <- function() {
    par(mfrow = c(1, 2), mar = c(5, 5.2, 4, 1), family = FONT,
        cex.main = 1.9, cex.lab = 1.8, cex.axis = 1.5, lwd = 1.2)
    for (ind in 1:2) {
      idx <- ((ind - 1) * 16 + 1):(ind * 16)
      matplot(1:16, t(mu[, idx]), type = "l", lty = 1, lwd = 3.5, col = cols,
              xlab = "time", ylab = "effect",
              main = sprintf("Individual %d", ind), ylim = range(mu))
      if (ind == 1) legend("topleft", legend = paste("cluster", 1:4),
                           col = cols, lwd = 3.5, bty = "n", cex = 1.5)
    }
  }
  save_base("fig_1B_mirror_mechanism", draw, w = 14, h = 6); message("[fig] 1B mirror mechanism done")
}

## ---- Fig 1D: trace plots (base graphics) ----
d_files <- list.files(file.path(MCG_OUT_DIR, "sim1D"), pattern = "^conv", full.names = TRUE)
if (length(d_files)) {
  x <- readRDS(d_files[1]); ml <- x$chains
  pars <- intersect(c("phi","phi_free","p[1]","p[2]","beta[1, 1]"), coda::varnames(ml))
  if (length(pars)) {
    draw <- function() {
      par(mfrow = c(length(pars), 1), mar = c(4.5, 5.2, 3, 1), family = FONT,
          cex.main = 1.8, cex.lab = 1.7, cex.axis = 1.4)
      for (pn in pars) {
        rng <- range(sapply(ml, function(c) c[, pn]))
        plot(NA, xlim = c(1, nrow(ml[[1]])), ylim = rng, xlab = "iter", ylab = pn,
             main = paste("trace:", pn))
        for (ci in seq_along(ml)) lines(ml[[ci]][, pn], col = ci, lwd = 1.4)
      }
    }
    save_base("fig_1D_trace", draw, w = 12, h = 11); message("[fig] 1D trace done")
  }
}
message("[03-local] figures in ", FIG)

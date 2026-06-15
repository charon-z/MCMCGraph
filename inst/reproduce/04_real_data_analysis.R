#!/usr/bin/env Rscript
# 04_real_data_analysis.R -------------------------------------------------
# End-to-end real-data workflow with the published MCMCGraph functions:
#   1. load a paired longitudinal matrix (n x 2d)
#   2. BIC-based selection of the number of clusters J  (fit_many_J / eval_bic)
#   3. final clustering at the selected J
#   4. cluster mean-curve table + posterior cluster probabilities
#   5. trace/density diagnostics for key parameters
#   6. a module-level network graph (clusters linked by mean-curve similarity)
#
# Dataset selection (environment variable MCG_DATASET):
#   "example" (default) : the bundled example_binary toy data (runs in minutes)
#   "lincs"             : LINCS L1000 vorinostat MCF7/PC3 (978 x 12), if present
#   "smillie"           : Smillie UC pseudobulk (if present)
# Real datasets are large; see data/*/README and the scripts/ directory for the
# full preprocessing. The example path is what the Quick Start in the README
# exercises, so a reviewer can reproduce the workflow without the big files.

suppressPackageStartupMessages(library(MCMCGraph))
src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "MCMCGraph/inst/reproduce"
source(file.path(src_dir, "config.R"))
OUT <- file.path(MCG_OUT_DIR, "real"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dataset <- Sys.getenv("MCG_DATASET", unset = "example")
J_grid  <- eval(parse(text = Sys.getenv("MCG_JGRID", unset = "2:6")))
niter   <- as.integer(Sys.getenv("MCG_NITER", unset = "8000"))

load_dataset <- function(name) {
  if (name == "lincs") {
    f <- file.path(MCG_ROOT, "data/LINCS_L1000_vorinostat/processed",
                   "LINCS_Vorinostat_MCF7_PC3_978x12_expression.tsv")
    if (!file.exists(f)) stop("LINCS file not found: ", f)
    df <- read.delim(f, check.names = FALSE, row.names = 1)
    df$gene_symbol <- NULL
    return(list(Y = as.matrix(df), times = 1:6, label = "LINCS vorinostat (MCF7|PC3)"))
  }
  if (name == "smillie") {
    f <- file.path(MCG_ROOT, "data/SmillieUC/processed", "SmillieUC_5000x20_expression.tsv")
    if (!file.exists(f)) stop("Smillie file not found: ", f)
    Y <- as.matrix(read.delim(f, check.names = FALSE, row.names = 1))
    return(list(Y = Y, times = seq_len(ncol(Y) / 2), label = "Smillie UC"))
  }
  data("example_binary", package = "MCMCGraph", envir = environment())
  list(Y = example_binary$y, times = example_binary$times, label = "example_binary (toy)")
}

ds <- load_dataset(dataset)
message(sprintf("[04] dataset=%s  n=%d  d=%d", ds$label, nrow(ds$Y), ncol(ds$Y)))

## 1-2. BIC sweep over J --------------------------------------------------
sweep <- MCMCGraph::fit_many_J(ds$Y, J_grid = J_grid, times = ds$times, niter = niter)
ev <- sweep$eval
write.csv(ev, file.path(OUT, sprintf("%s_BIC.csv", dataset)), row.names = FALSE)
bestJ <- ev$J[which.min(ev$BIC)]
message(sprintf("[04] BIC-selected J = %d", bestJ))

pdf(file.path(OUT, sprintf("%s_BIC.pdf", dataset)), width = 6, height = 4.5)
MCMCGraph::plot_bic(ev); dev.off()

## 3. final clustering at best J ------------------------------------------
fit <- sweep$fits[[paste0("J", bestJ)]]
if (is.null(fit)) fit <- MCMCGraph::run_mcmc_binary(ds$Y, J = bestJ, times = ds$times, niter = niter)

clusters <- data.frame(feature = rownames(ds$Y), cluster = fit$clustering,
                       fit$cluster_prob, uncertainty = fit$cluster_uncertainty)
write.csv(clusters, file.path(OUT, sprintf("%s_J%d_clusters.csv", dataset, bestJ)), row.names = FALSE)

## 4. cluster mean curves -------------------------------------------------
d_single <- fit$model_info$d_single
Z0 <- MCMCGraph:::make_Z0_binary(ds$times)
mu <- fit$posterior_mean$beta %*% t(Z0)        # J x 2d
write.csv(data.frame(cluster = seq_len(bestJ), mu),
          file.path(OUT, sprintf("%s_J%d_mean_curves.csv", dataset, bestJ)), row.names = FALSE)

## 5. trace/density diagnostics -------------------------------------------
pdf(file.path(OUT, sprintf("%s_J%d_trace.pdf", dataset, bestJ)), width = 7, height = 7)
MCMCGraph::plot_trace_density(fit, params = c("^phi", "^p\\[")); dev.off()

## 6. module-level network graph (clusters linked by mean-curve similarity)-
S <- cor(t(mu))                                 # J x J correlation of mean curves
adj <- S; diag(adj) <- 0; adj[abs(adj) < 0.5] <- 0
pdf(file.path(OUT, sprintf("%s_J%d_module_network.pdf", dataset, bestJ)), width = 6, height = 6)
if (requireNamespace("igraph", quietly = TRUE)) {
  g <- igraph::graph_from_adjacency_matrix(abs(adj), mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::V(g)$size <- 8 + 30 * as.numeric(table(factor(fit$clustering, levels = 1:bestJ))) / nrow(ds$Y)
  ecol <- ifelse(adj[igraph::as_edgelist(g, names = FALSE)] > 0, "#d95f02", "#1b9e77")
  plot(g, edge.width = 2 * igraph::E(g)$weight, edge.color = ecol,
       vertex.label = paste0("M", 1:bestJ), main = "Module similarity network")
} else {
  image(seq_len(bestJ), seq_len(bestJ), S, xlab = "module", ylab = "module",
        main = "Module mean-curve correlation")
}
dev.off()
message("[04] outputs written to ", OUT)

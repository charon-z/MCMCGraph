#!/usr/bin/env Rscript

env_value <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

method <- tolower(env_value("BPFC_BENCH_METHOD", "bpfc"))
dataset <- tolower(env_value("BPFC_BENCH_DATASET", "sim"))
n_requested <- as.integer(env_value("BPFC_BENCH_N", "3000"))
J <- as.integer(env_value("BPFC_BENCH_J", "8"))
niter <- as.integer(env_value("BPFC_BENCH_NITER", "5000"))
seed <- as.integer(env_value("BPFC_BENCH_SEED", "20260904"))
repeat_id <- as.integer(env_value("BPFC_BENCH_REPEAT", "1"))
output_file <- env_value("BPFC_BENCH_OUTPUT", "benchmark_result.csv")
mouse_input <- env_value("BPFC_BENCH_MOUSE_INPUT", "")
source_dir <- env_value("BPFC_BENCH_SOURCE", "")

allowed_methods <- c("bpfc", "bpfc_em", "gmm", "kmeans", "agglomerative")
if (!method %in% allowed_methods) {
  stop("Unknown method: ", method)
}
if (!dataset %in% c("sim", "mouse")) {
  stop("BPFC_BENCH_DATASET must be 'sim' or 'mouse'.")
}
if (!is.finite(J) || J < 2L) stop("BPFC_BENCH_J must be at least 2.")
if (method == "bpfc" && (!is.finite(niter) || niter < 20L)) {
  stop("BPFC_BENCH_NITER must be at least 20 for BPFC-MCMC.")
}

suppressPackageStartupMessages(library(BPFC))

legendre_basis <- function(times, order = 4L) {
  scaled <- -1 + 2 * (times - min(times)) / (max(times) - min(times))
  full <- cbind(
    1,
    scaled,
    0.5 * (3 * scaled^2 - 1),
    0.5 * (5 * scaled^3 - 3 * scaled),
    0.125 * (35 * scaled^4 - 30 * scaled^2 + 3)
  )
  full[, seq_len(order + 1L), drop = FALSE]
}

joint_basis <- function(d_single = 16L) {
  one <- legendre_basis(seq_len(d_single))
  q <- ncol(one)
  basis <- matrix(0, nrow = 2L * d_single, ncol = 2L * q)
  basis[seq_len(d_single), seq_len(q)] <- one
  basis[d_single + seq_len(d_single), q + seq_len(q)] <- one
  basis
}

generate_simulated_data <- function(n, J, seed, d_single = 16L) {
  set.seed(seed)
  basis <- joint_basis(d_single)
  coefficient_index <- seq_len(J * ncol(basis))
  beta <- matrix(
    0.55 * sin(0.43 * coefficient_index) +
      0.25 * cos(0.19 * coefficient_index),
    nrow = J, byrow = TRUE
  )
  beta[, c(1L, 6L)] <- beta[, c(1L, 6L)] +
    seq(-1.1, 1.1, length.out = J)
  component_mean <- beta %*% t(basis)

  labels <- rep(seq_len(J), length.out = n)
  labels <- sample(labels, length(labels), replace = FALSE)
  innovations <- matrix(stats::rnorm(n * 2L * d_single, sd = 0.45), nrow = n)
  residual <- innovations
  phi <- 0.55
  for (block_start in c(1L, d_single + 1L)) {
    block_end <- block_start + d_single - 1L
    if (block_end > block_start) {
      for (time_index in (block_start + 1L):block_end) {
        residual[, time_index] <- innovations[, time_index] -
          phi * innovations[, time_index - 1L]
      }
    }
  }
  y <- component_mean[labels, , drop = FALSE] + residual
  colnames(y) <- c(
    sprintf("state1_t%02d", seq_len(d_single)),
    sprintf("state2_t%02d", seq_len(d_single))
  )
  list(y = y, truth = labels, times = seq_len(d_single))
}

lop_features <- function(y, times) {
  basis <- joint_basis(length(times))
  inverse_crossproduct <- solve(crossprod(basis) + diag(1e-8, ncol(basis)))
  t(inverse_crossproduct %*% crossprod(basis, t(y)))
}

if (dataset == "sim") {
  if (!is.finite(n_requested) || n_requested < J) {
    stop("BPFC_BENCH_N must be at least J for simulated data.")
  }
  data_object <- generate_simulated_data(n_requested, J, seed)
  y <- data_object$y
  truth <- data_object$truth
  times <- data_object$times
} else {
  if (!nzchar(mouse_input) || !file.exists(mouse_input)) {
    stop("A valid BPFC_BENCH_MOUSE_INPUT file is required for mouse data.")
  }
  y <- as.matrix(utils::read.table(mouse_input, header = TRUE, check.names = FALSE))
  storage.mode(y) <- "double"
  if (ncol(y) %% 2L != 0L) stop("Mouse input must contain two equal time blocks.")
  truth <- NULL
  times <- seq_len(ncol(y) / 2L)
}

n <- nrow(y)
d <- ncol(y)
gc(reset = TRUE)
start_time <- proc.time()[["elapsed"]]

fit <- switch(
  method,
  bpfc = BPFC::run_mcmc_binary(
    y, J = J, format = "wide", times = times,
    niter = niter, thin = 1L, burnin_frac = 0.25,
    seed = seed
  ),
  bpfc_em = {
    if (!nzchar(source_dir)) {
      stop("BPFC_BENCH_SOURCE is required for the BPFC-EM comparator.")
    }
    helper <- file.path(source_dir, "inst", "reproduce", "lib_sim.R")
    if (!file.exists(helper)) stop("Cannot find BPFC-EM helper: ", helper)
    source(helper, local = .GlobalEnv)
    basis <- mcg_basis_joint(length(times), order = 4L, times = times)
    mcg_em_sad(
      y, J = J,
      blocks = list(seq_along(times), length(times) + seq_along(times)),
      basis = basis, n_init = 10L, seed = seed,
      max_iter = 200L, tol = 1e-5
    )
  },
  gmm = {
    if (!requireNamespace("mclust", quietly = TRUE)) {
      stop("Package 'mclust' is required for the GMM comparator.")
    }
    suppressPackageStartupMessages(library(mclust))
    features <- scale(lop_features(y, times))
    mclust::Mclust(features, G = J, modelNames = "VVV", verbose = FALSE)
  },
  kmeans = {
    set.seed(seed)
    features <- scale(lop_features(y, times))
    stats::kmeans(features, centers = J, nstart = 100L, iter.max = 100L)
  },
  agglomerative = {
    features <- scale(lop_features(y, times))
    tree <- stats::hclust(stats::dist(features), method = "ward.D2")
    list(cluster = stats::cutree(tree, k = J), tree = tree)
  }
)

elapsed_seconds <- proc.time()[["elapsed"]] - start_time
predicted <- switch(
  method,
  bpfc = fit$clustering,
  bpfc_em = fit$z,
  gmm = fit$classification,
  kmeans = fit$cluster,
  agglomerative = fit$cluster
)

ari <- NA_real_
if (!is.null(truth) && requireNamespace("mclust", quietly = TRUE)) {
  ari <- mclust::adjustedRandIndex(truth, predicted)
}

result <- data.frame(
  run_id = tools::file_path_sans_ext(basename(output_file)),
  method = method,
  dataset = dataset,
  n = n,
  dimensions = d,
  J = J,
  niter = if (method == "bpfc") niter else NA_integer_,
  repeat_id = repeat_id,
  seed = seed,
  elapsed_seconds = unname(elapsed_seconds),
  result_object_mb = as.numeric(utils::object.size(fit)) / 1024^2,
  adjusted_rand_index = ari,
  bpfc_version = as.character(utils::packageVersion("BPFC")),
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  stringsAsFactors = FALSE
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, output_file, row.names = FALSE)
print(result)

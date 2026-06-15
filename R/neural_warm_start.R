#' Build a neural-network warm start for paired Bayesian MCMC clustering
#'
#' This helper trains a small single-hidden-layer neural classifier on labelled
#' simulated paired trajectories, predicts cluster probabilities for target
#' trajectories, and converts the predicted labels into the initialization list
#' expected by [run_mcmc_binary()]. The neural network only chooses the starting
#' point; posterior inference is still performed by the Bayesian MCMC model.
#'
#' @param train_x labelled training trajectories in the same format as
#'   [as_binary_data()].
#' @param train_z integer training labels in `1:J`.
#' @param target_x target trajectories to initialize.
#' @param J number of clusters. Defaults to `max(train_z)`.
#' @param train_format,target_format data formats passed to [as_binary_data()].
#' @param times optional time grid for each individual.
#' @param hidden number of hidden units in the neural initializer.
#' @param decay weight decay passed to `nnet::nnet`.
#' @param maxit maximum neural-network optimization iterations.
#' @param seed random seed for neural-network training.
#' @param trace whether `nnet` should print training progress.
#' @return An initialization list with `z`, `beta`, `v_sq`, `phi` and `p`, plus
#'   neural predicted probabilities and training metadata.
#' @export
neural_warmstart_binary <- function(
    train_x,
    train_z,
    target_x,
    J = NULL,
    train_format = c("auto", "wide", "long"),
    target_format = c("auto", "wide", "long"),
    times = NULL,
    hidden = 12L,
    decay = 1e-4,
    maxit = 300L,
    seed = 123,
    trace = FALSE
) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Package 'nnet' is required for neural warm starts.", call. = FALSE)
  }
  train_format <- match.arg(train_format)
  target_format <- match.arg(target_format)
  set.seed(seed)

  train <- as_binary_data(train_x, format = train_format, times = times)
  target <- as_binary_data(target_x, format = target_format, times = train$times_single)
  if (train$d_single != target$d_single) {
    stop("train_x and target_x must have the same number of time points.", call. = FALSE)
  }

  train_z <- as.integer(train_z)
  if (length(train_z) != nrow(train$y)) {
    stop("train_z must have length nrow(train_x).", call. = FALSE)
  }
  if (is.null(J)) J <- max(train_z)
  J <- as.integer(J)
  if (any(train_z < 1L | train_z > J)) {
    stop("train_z must contain integer labels in 1:J.", call. = FALSE)
  }

  x_train <- paired_neural_features_binary(train$y, train$times_single)
  x_target <- paired_neural_features_binary(target$y, target$times_single)
  center <- colMeans(x_train)
  scale <- apply(x_train, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  x_train <- scale(x_train, center = center, scale = scale)
  x_target <- scale(x_target, center = center, scale = scale)

  y_train <- nnet::class.ind(factor(train_z, levels = seq_len(J)))
  max_wts <- (ncol(x_train) + 1L) * hidden + (hidden + 1L) * J + 100L
  fit <- nnet::nnet(
    x = x_train,
    y = y_train,
    size = hidden,
    softmax = TRUE,
    decay = decay,
    maxit = maxit,
    trace = trace,
    MaxNWts = max_wts
  )
  prob <- stats::predict(fit, x_target, type = "raw")
  prob <- as.matrix(prob)
  if (ncol(prob) != J) {
    stop("Neural prediction returned an unexpected number of columns.", call. = FALSE)
  }
  colnames(prob) <- paste0("cluster", seq_len(J))
  z0 <- max.col(prob, ties.method = "first")

  Z0 <- make_Z0_binary(target$times_single)
  init <- initialization_from_labels_binary(target$y, z0, J, Z0)
  init$method <- "neural_warmstart"
  init$neural_prob <- prob
  init$neural_confidence <- apply(prob, 1L, max)
  init$neural_info <- list(
    hidden = hidden,
    decay = decay,
    maxit = maxit,
    seed = seed,
    train_n = nrow(train$y),
    target_n = nrow(target$y),
    n_features = ncol(x_train)
  )
  init
}

#' Run Bayesian MCMC with a neural-network warm start
#'
#' This is a convenience wrapper around `neural_warmstart_binary()` and
#' [run_mcmc_binary()]. It leaves the original MCMC model unchanged; the neural
#' network supplies only the initial labels and parameter values.
#'
#' @param x target trajectories.
#' @param J number of clusters.
#' @param train_x labelled simulated training trajectories.
#' @param train_z integer simulated training labels.
#' @param ... arguments passed to [run_mcmc_binary()], such as `niter`, `thin`,
#'   `n_batch`, `burnin_frac`, `priors` and `seed`.
#' @param format target data format passed to [run_mcmc_binary()].
#' @param times optional target time grid.
#' @param train_format training data format.
#' @param hidden,decay,maxit,trace neural initializer controls.
#' @param neural_seed random seed for neural-network training.
#' @return An `mcmcgraph_result` with an added `warmstart` element.
#' @export
run_mcmc_binary_neural_warmstart <- function(
    x,
    J,
    train_x,
    train_z,
    ...,
    format = c("auto", "wide", "long"),
    times = NULL,
    train_format = c("auto", "wide", "long"),
    hidden = 12L,
    decay = 1e-4,
    maxit = 300L,
    neural_seed = 123,
    trace = FALSE
) {
  format <- match.arg(format)
  train_format <- match.arg(train_format)
  init <- neural_warmstart_binary(
    train_x = train_x,
    train_z = train_z,
    target_x = x,
    J = J,
    train_format = train_format,
    target_format = format,
    times = times,
    hidden = hidden,
    decay = decay,
    maxit = maxit,
    seed = neural_seed,
    trace = trace
  )
  fit <- run_mcmc_binary(
    x = x,
    J = J,
    format = format,
    times = times,
    init = init,
    ...
  )
  fit$warmstart <- list(
    method = init$method,
    z = init$z,
    prob = init$neural_prob,
    confidence = init$neural_confidence,
    info = init$neural_info
  )
  fit
}

#' @keywords internal
paired_neural_features_binary <- function(y, times_single = NULL) {
  y <- as.matrix(y)
  if (!is.numeric(y)) storage.mode(y) <- "double"
  if (ncol(y) %% 2L != 0L) stop("y must have 2*d columns.", call. = FALSE)
  d_single <- ncol(y) / 2L
  if (is.null(times_single)) times_single <- seq_len(d_single)
  Z0 <- make_Z0_binary(times_single)
  beta <- t(solve(crossprod(Z0) + diag(1e-8, ncol(Z0)), crossprod(Z0, t(y))))
  q <- ncol(beta) / 2L
  b1 <- beta[, seq_len(q), drop = FALSE]
  b2 <- beta[, q + seq_len(q), drop = FALSE]
  raw_diff <- y[, seq_len(d_single), drop = FALSE] -
    y[, d_single + seq_len(d_single), drop = FALSE]
  cbind(
    y,
    raw_diff,
    beta,
    b1 - b2,
    abs(b1 - b2),
    b1 * b2
  )
}

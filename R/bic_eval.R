#' Fit multiple J values and return BIC table
#'
#' @param x input data (wide/long)
#' @param J_grid integer vector of candidate J
#' @param ... passed to run_mcmc_binary (e.g., niter, thin)
#' @return list(fits, eval) where eval is data.frame(J,BIC,loglik)
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' res <- fit_many_J(example_binary$y, J_grid = 2:4,
#'                   times = example_binary$times, niter = 300)
#' res$eval
#' }
fit_many_J <- function(x, J_grid = 2:6, ...) {
  fits <- vector("list", length(J_grid))
  evals <- vector("list", length(J_grid))

  for (i in seq_along(J_grid)) {
    J <- J_grid[i]
    fit <- run_mcmc_binary(x, J = J, ...)
    ev <- eval_bic(fit, x = x, format = fit$data_info$format %||% "auto")
    fits[[i]] <- fit
    evals[[i]] <- ev
  }
  names(fits) <- paste0("J", J_grid)
  list(fits = fits, eval = do.call(rbind, evals))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Compute BIC from a fitted object (posterior mean plug-in)
#'
#' Uses the relabeled posterior means stored on the fit (so the plug-in is not
#' corrupted by label switching). The two individuals share a single SAD
#' time-correlation parameter phi (paper 2.3).
#'
#' @param fit mcmcgraph_result
#' @param x original input data (same used in fitting)
#' @param format "auto" | "wide" | "long"
#' @return data.frame with J, loglik, BIC
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' fit <- run_mcmc_binary(example_binary$y, J = 3, times = example_binary$times,
#'                        niter = 300, seed = 1)
#' eval_bic(fit, x = example_binary$y, format = "wide")
#' }
eval_bic <- function(fit, x, format = c("auto","wide","long")) {
  format <- match.arg(format)
  dat <- as_binary_data(x, format = format, times = fit$data_info$times_single %||% NULL)
  y <- dat$y
  n <- nrow(y)
  d <- ncol(y)
  d_single <- dat$d_single

  J <- fit$model_info$J
  P <- fit$model_info$P
  Z0_binary <- make_Z0_binary(dat$times_single)

  # relabeled posterior-mean plug-in estimates
  pm <- fit$posterior_mean
  phi_hat  <- as.numeric(pm$phi)
  p_hat    <- pmax(as.numeric(pm$p), 1e-12); p_hat <- p_hat / sum(p_hat)
  v_hat    <- pmax(as.numeric(pm$v_sq), 1e-8)
  beta_hat <- pm$beta                       # J x P

  mu_mat <- beta_hat %*% t(Z0_binary)       # J x d

  # log density for the 2-block SAD recursion with a single shared phi
  dSAD_log_2block <- function(xi, mean, v_sq, phi, d_single) {
    q <- 0
    s <- 0
    for (t in 1:d_single) {
      s <- (xi[t] - mean[t]) + phi * s
      q <- q + (s * s) / v_sq[t] + log(2 * pi * v_sq[t])
    }
    s <- 0
    for (t in (d_single + 1):(2 * d_single)) {
      s <- (xi[t] - mean[t]) + phi * s
      q <- q + (s * s) / v_sq[t] + log(2 * pi * v_sq[t])
    }
    -0.5 * q
  }

  loglik <- 0
  for (i in 1:n) {
    yi <- as.numeric(y[i, ])
    log_comp <- numeric(J)
    for (j in 1:J) {
      logf <- dSAD_log_2block(yi, mu_mat[j, ], v_hat, phi_hat, d_single)
      log_comp[j] <- log(p_hat[j]) + logf
    }
    loglik <- loglik + log_sum_exp(log_comp)
  }

  # free parameters (fixed order = 4, K = 2):
  # beta: J*P, v_sq: 2*d_single, phi: 1 (shared), mixing weights: J-1
  k <- J * P + d + 1 + (J - 1)
  bic <- -2 * loglik + k * log(n)
  data.frame(J = J, loglik = loglik, BIC = bic)
}

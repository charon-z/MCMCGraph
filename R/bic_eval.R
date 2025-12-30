#' Fit multiple J values and return BIC table
#'
#' @param x input data (wide/long)
#' @param J_grid integer vector of candidate J
#' @param ... passed to run_mcmc_binary (e.g., niter, thin)
#' @return list(fits, eval) where eval is data.frame(J,BIC,loglik)
#' @export
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
#' @param fit mcmcgraph_result
#' @param x original input data (same used in fitting)
#' @param format "auto" | "wide" | "long"
#' @return data.frame with J, loglik, BIC
#' @export
eval_bic <- function(fit, x, format = c("auto","wide","long")) {
  format <- match.arg(format)
  dat <- as_binary_data(x, format = format, times = fit$data_info$times_single %||% NULL)
  y <- dat$y
  n <- nrow(y)
  d <- ncol(y)
  d_single <- dat$d_single

  # posterior mean plug-in
  samp <- fit$posterior_samples$all_params
  mean_vec <- colMeans(samp)
  nm <- names(mean_vec)
  J <- fit$model_info$J
  P <- fit$model_info$P
  Z0_binary <- make_Z0_binary(dat$times_single)

  phi1_hat <- as.numeric(mean_vec["phi1"])
  phi2_hat <- as.numeric(mean_vec["phi2"])

  # p
  p_names <- grep("^p\\[", nm, value = TRUE)
  p_hat <- numeric(J)
  jidx <- as.integer(sub("^p\\[(\\d+)\\]$", "\\1", p_names))
  p_hat[jidx] <- as.numeric(mean_vec[p_names])
  p_hat <- pmax(p_hat, 1e-12); p_hat <- p_hat / sum(p_hat)

  # v_sq
  v_names <- grep("^v_sq\\[", nm, value = TRUE)
  v_hat <- numeric(d)
  tidx <- as.integer(sub("^v_sq\\[(\\d+)\\]$", "\\1", v_names))
  v_hat[tidx] <- as.numeric(mean_vec[v_names])
  v_hat <- pmax(v_hat, 1e-8)

  # beta matrix J x P
  b_names <- grep("^beta\\[", nm, value = TRUE)
  beta_hat <- matrix(0, nrow = J, ncol = P)
  rx <- regexec("^beta\\[(\\d+),\\s*(\\d+)\\]$", b_names)
  parsed <- regmatches(b_names, rx)
  for (ii in seq_along(b_names)) {
    jj <- as.integer(parsed[[ii]][2])
    kk <- as.integer(parsed[[ii]][3])
    beta_hat[jj, kk] <- as.numeric(mean_vec[b_names[ii]])
  }

  mu_mat <- beta_hat %*% t(Z0_binary)  # J x d

  # log density for 2-block SAD recursion
  dSAD_log_2block <- function(xi, mean, v_sq, phi1, phi2, d_single) {
    q <- 0
    s <- 0
    for (t in 1:d_single) {
      s <- (xi[t] - mean[t]) + phi1 * s
      q <- q + (s * s) / v_sq[t] + log(2 * pi * v_sq[t])
    }
    s <- 0
    for (t in (d_single + 1):(2 * d_single)) {
      s <- (xi[t] - mean[t]) + phi2 * s
      q <- q + (s * s) / v_sq[t] + log(2 * pi * v_sq[t])
    }
    -0.5 * q
  }

  loglik <- 0
  for (i in 1:n) {
    yi <- as.numeric(y[i, ])
    log_comp <- numeric(J)
    for (j in 1:J) {
      logf <- dSAD_log_2block(yi, mu_mat[j, ], v_hat, phi1_hat, phi2_hat, d_single)
      log_comp[j] <- log(p_hat[j]) + logf
    }
    loglik <- loglik + log_sum_exp(log_comp)
  }

  # param count for your fixed order=4, K=2 setup:
  # beta: J*10, v_sq: 2d_single, phi:2, p:(J-1)
  k <- J * P + d + 2 + (J - 1)
  bic <- -2 * loglik + k * log(n)
  data.frame(J = J, loglik = loglik, BIC = bic)
}

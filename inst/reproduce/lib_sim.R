# lib_sim.R ---------------------------------------------------------------
# Shared simulation / inference machinery for the BPFC reproduce suite.
#
# Provides:
#   * Legendre (LOP) basis at arbitrary order and block layout.
#   * A model-consistent two-block SAD(1) data generator (single shared phi).
#   * Faithful generators for the three paper scenarios (K = 3 / 5 / 8), reading
#     the stored truth parameters in benchmark/data/sim_truth_K*.rds.
#   * A "mirror" scenario in which clusters are distinguishable only through the
#     pairing of the two individuals (the litmus test for joint clustering).
#   * A from-scratch EM for the SAD(1)+LOP finite mixture (the *same* model as
#     the MCMC), with multi-start initialisation -> used for the EM-vs-MCMC
#     stability and init-sensitivity comparisons.
#   * A generalised nimble MCMC runner over an arbitrary block structure and
#     basis (so joint / concatenation / independent share one kernel and only
#     the "clustering unit" changes -> a clean ablation).
#   * Baseline method wrappers (LOP k-means, SVM-RBF, XGBoost).
#
# Depends: nimble, BPFC (for .relabel_ecr), mclust, clue; e1071/xgboost
# only when the supervised wrappers are called.

suppressPackageStartupMessages({
  library(nimble)
  library(BPFC)
})

# ---- Legendre basis -------------------------------------------------------

# Shifted Legendre polynomials up to `order` evaluated at `times` -> d x (order+1).
mcg_legendre <- function(times, order = 4L) {
  d <- length(times)
  tmin <- min(times); tmax <- max(times)
  tx <- if (tmax == tmin) rep(0, d) else (-1 + 2 * (times - tmin) / (tmax - tmin))
  if (order > 5L) stop("mcg_legendre supports order <= 5")
  P6 <- function(x) c(1, x, 0.5 * (3 * x^2 - 1), 0.5 * (5 * x^3 - 3 * x),
                      0.125 * (35 * x^4 - 30 * x^2 + 3),
                      0.125 * (63 * x^5 - 70 * x^3 + 15 * x))
  Z <- t(vapply(tx, P6, numeric(6)))
  Z[, seq_len(order + 1L), drop = FALSE]
}

# Block-diagonal "joint" paired basis: 2d x 2(order+1).
mcg_basis_joint <- function(d_single, order = 4L, times = seq_len(d_single)) {
  Z1 <- mcg_legendre(times, order); q <- ncol(Z1)
  Zb <- matrix(0, 2 * d_single, 2 * q)
  Zb[1:d_single, 1:q] <- Z1
  Zb[(d_single + 1L):(2 * d_single), (q + 1L):(2 * q)] <- Z1
  Zb
}

# Single continuous basis over `d` points: d x (order+1).
mcg_basis_single <- function(d, order = 4L, times = seq_len(d)) {
  mcg_legendre(times, order)
}

# ---- model-consistent SAD generator --------------------------------------

# Generate one paired vector consistent with dSADmvnorm: within each block
# x_t = m_t + e_t - phi * e_{t-1}, e_t ~ N(0, v_t), e_0 = 0. Shared phi.
mcg_gen_two_block <- function(mu, phi, v_sq, d_single) {
  d <- 2L * d_single
  out <- numeric(d)
  for (blk in 0:1) {
    e_prev <- 0
    for (t in 1:d_single) {
      idx <- blk * d_single + t
      e <- rnorm(1, 0, sqrt(v_sq[idx]))
      out[idx] <- mu[idx] + e - phi * e_prev
      e_prev <- e
    }
  }
  out
}

# ---- paper scenarios ------------------------------------------------------

# Locate benchmark/data relative to a base directory.
mcg_truth_path <- function(K, data_dir) file.path(data_dir, sprintf("sim_truth_K%d.rds", K))

# Return the stored generative parameters for a scenario (K in {3,5,8}).
mcg_scenario_params <- function(K, data_dir) {
  s <- readRDS(mcg_truth_path(K, data_dir))
  list(
    K = as.integer(s$K_true),
    n = as.integer(s$n),
    d_single = 16L,
    times = if (!is.null(s$Times)) as.numeric(s$Times) else 1:16,
    p = as.numeric(s$p_true),
    phi = as.numeric(s$phi1_true),          # single shared phi (paper method)
    v_sq = as.numeric(s$v_sq_true),
    beta = s$beta_true                        # K x 10
  )
}

# Generate a fresh replicate of a paper scenario under the single-phi model.
mcg_gen_scenario <- function(par, seed) {
  set.seed(seed)
  K <- par$K; n <- par$n; d_single <- par$d_single
  Z0 <- mcg_basis_joint(d_single, order = 4L, times = par$times)
  mu <- par$beta %*% t(Z0)                    # K x 2d
  z <- sample.int(K, n, replace = TRUE, prob = par$p)
  Y <- matrix(0, n, 2L * d_single)
  for (i in 1:n) Y[i, ] <- mcg_gen_two_block(mu[z[i], ], par$phi, par$v_sq, d_single)
  colnames(Y) <- c(sprintf("i1_t%02d", 1:d_single), sprintf("i2_t%02d", 1:d_single))
  list(Y = Y, z_true = z, K = K, d_single = d_single, par = par)
}

# ---- mirror scenario ------------------------------------------------------
# The litmus test for joint clustering. Four clusters sit at mu0 +/- d1 +/- d2,
# where the displacement directions d1, d2 are chosen to lie in the *joint*
# (block-diagonal LOP) mean subspace but to be ORTHOGONAL to the
# *concatenation* subspace (a single LOP over the 2d-long concatenated axis).
#
# Consequence, by construction:
#   * Concatenation projects all four cluster means onto the SAME point -> it is
#     representationally blind to the differences and can only cluster on noise.
#   * The joint model keeps the two individuals on separate bases and sees the
#     full mu0 +/- d1 +/- d2 separation.
#   * A naive "cluster each individual into K then match modules 1-1"
#     (independent) approach also breaks, because each individual marginally
#     shows fewer groups than the K joint modules, so the 1-1 module matching is
#     ill-posed (the "active in ind1 / suppressed in ind2" motif).
# This encodes exactly the empirical Cluster 5/7 phenomenon: a module whose two
# individuals move in opposite directions cannot be recovered by collapsing or
# by per-individual matching, only by treating the pair as one unit.
mcg_gen_mirror <- function(seed, n = 2000L, d_single = 16L, phi = 0.6,
                           v_scale = 1.0, amp = 1.0, p = c(0.30, 0.30, 0.20, 0.20)) {
  set.seed(seed)
  times <- seq_len(d_single)
  Zj <- mcg_basis_joint(d_single, 4L, times)     # 2d x 10  (joint / block-diagonal)
  Zc <- mcg_basis_single(2L * d_single, 4L)       # 2d x 5   (concatenation basis)
  # directions in col(Zj) that are orthogonal to col(Zc):
  #   delta = Zj %*% c  with  t(Zc) %*% delta = 0  ->  (t(Zc) %*% Zj) c = 0.
  M <- crossprod(Zc, Zj)                           # 5 x 10
  ns <- svd(M, nu = 0, nv = ncol(M))
  sv_full <- c(ns$d, rep(0, ncol(M) - length(ns$d)))   # pad singular values to 10
  null_dirs <- ns$v[, which(sv_full < 1e-8), drop = FALSE]
  stopifnot(ncol(null_dirs) >= 2)
  # pick two null directions whose induced mean displacements are well separated
  d1 <- as.numeric(Zj %*% null_dirs[, 1]); d2 <- as.numeric(Zj %*% null_dirs[, 2])
  d1 <- amp * d1 / sqrt(mean(d1^2)); d2 <- amp * d2 / sqrt(mean(d2^2))
  beta0 <- c(0.5, 0.6, 0.2, 0, 0, 0.5, 0.6, 0.2, 0, 0)
  mu0 <- as.numeric(Zj %*% beta0)
  signs <- list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  K <- length(signs)
  mu <- t(vapply(signs, function(s) mu0 + s[1] * d1 + s[2] * d2, numeric(2 * d_single)))
  v_sq <- rep(v_scale * 0.25, 2 * d_single)
  z <- sample.int(K, n, replace = TRUE, prob = p)
  Y <- matrix(0, n, 2 * d_single)
  for (i in 1:n) Y[i, ] <- mcg_gen_two_block(mu[z[i], ], phi, v_sq, d_single)
  colnames(Y) <- c(sprintf("i1_t%02d", 1:d_single), sprintf("i2_t%02d", 1:d_single))
  list(Y = Y, z_true = z, K = K, d_single = d_single,
       mu = mu, d1 = d1, d2 = d2, phi = phi, v_sq = v_sq, p = p)
}

# Factorial "pairing" mirror: the empirical Cluster 5/7 motif. Two individual
# shapes A (rising) and B (falling); the four clusters are the four pairings
# (A,A) (A,B) (B,A) (B,B). Each individual marginally shows only TWO shapes, so
# clustering one individual and matching modules 1-1 cannot represent a feature
# that is "shape A in individual 1 but shape B in individual 2" -> the naive
# independent pipeline must merge clusters. The joint unit keeps the pairing.
mcg_gen_mirror_pairing <- function(seed, n = 2000L, d_single = 16L, phi = 0.6,
                                   v_scale = 1.0, p = c(0.30, 0.30, 0.20, 0.20)) {
  set.seed(seed)
  times <- seq_len(d_single)
  shapeA <- c(0.5,  0.9, 0.2, 0, 0)   # rising
  shapeB <- c(0.5, -0.8, 0.3, 0, 0)   # falling
  combos <- list(c("A","A"), c("A","B"), c("B","A"), c("B","B"))
  shp <- list(A = shapeA, B = shapeB); K <- length(combos)
  beta <- t(vapply(combos, function(cb) c(shp[[cb[1]]], shp[[cb[2]]]), numeric(10)))
  Z0 <- mcg_basis_joint(d_single, 4L, times); mu <- beta %*% t(Z0)
  v_sq <- rep(v_scale * 0.12, 2 * d_single)
  z <- sample.int(K, n, replace = TRUE, prob = p)
  Y <- matrix(0, n, 2 * d_single)
  for (i in 1:n) Y[i, ] <- mcg_gen_two_block(mu[z[i], ], phi, v_sq, d_single)
  colnames(Y) <- c(sprintf("i1_t%02d", 1:d_single), sprintf("i2_t%02d", 1:d_single))
  list(Y = Y, z_true = z, K = K, d_single = d_single, mu = mu, phi = phi, v_sq = v_sq, p = p)
}

# ---- LOP feature extraction (for k-means / EM init) -----------------------

mcg_lop_features <- function(Y, d_single = ncol(Y) / 2L, order = 4L,
                             times = seq_len(d_single)) {
  Z0 <- mcg_basis_joint(d_single, order, times)
  XtXinv <- solve(crossprod(Z0) + diag(1e-8, ncol(Z0)))
  beta <- t(XtXinv %*% crossprod(Z0, t(Y)))
  colnames(beta) <- c(paste0("i1_b", seq_len(order + 1L)),
                      paste0("i2_b", seq_len(order + 1L)))
  beta
}

# =========================================================================
# EM for the SAD(1) + LOP finite mixture (same model as the MCMC).
# =========================================================================
# blocks: list of column-index vectors; the SAD recursion resets at each block
# start and all blocks share a single phi. basis: design matrix (d x P) mapping
# beta_j (length P) to the component mean. v_mode: "pertime" (heteroscedastic,
# the SAD default) or "single" (homoscedastic AR(1)). fix_phi: if non-NULL,
# phi is held fixed (phi = 0 gives the independence / iid-innovation model).

# Whitened innovations s (n x d) for residual matrix r given phi and block resets.
.mcg_whiten <- function(r, phi, resets) {
  d <- ncol(r); s <- matrix(0, nrow(r), d)
  for (t in seq_len(d)) {
    if (resets[t]) s[, t] <- r[, t] else s[, t] <- r[, t] + phi * s[, t - 1L]
  }
  s
}

# log SAD density of each row of Y under one component (mean vector mu).
.mcg_sad_loglik_rows <- function(Y, mu, phi, v_sq, resets) {
  r <- sweep(Y, 2, mu, "-")
  s <- .mcg_whiten(r, phi, resets)
  -0.5 * (rowSums(sweep(s^2, 2, v_sq, "/")) + sum(log(2 * pi * v_sq)))
}

mcg_em_sad <- function(Y, J, blocks, basis, order = 4L,
                       v_mode = c("pertime", "single"), fix_phi = NULL,
                       n_init = 10L, seed = 1L, max_iter = 200L, tol = 1e-5,
                       verbose = FALSE) {
  v_mode <- match.arg(v_mode)
  n <- nrow(Y); d <- ncol(Y); P <- ncol(basis)
  resets <- logical(d)
  for (b in blocks) resets[b[1]] <- TRUE
  ZtZ <- crossprod(basis)

  # GLS solve for beta given phi, v: minimise weighted whitened residual.
  # Build whitening operator columns of basis too, then weighted LS.
  solve_beta <- function(ybar, phi, v_sq) {
    # whiten basis columns and ybar with the same recursion, weight by 1/sqrt(v)
    Bw <- .mcg_whiten(t(basis), phi, resets)      # P x d  -> treat rows
    Bw <- t(Bw)                                   # d x P
    yw <- as.numeric(.mcg_whiten(matrix(ybar, 1), phi, resets))
    W <- 1 / v_sq
    lhs <- t(Bw * W) %*% Bw + diag(1e-8, P)
    rhs <- t(Bw * W) %*% yw
    as.numeric(solve(lhs, rhs))
  }

  one_run <- function(init_seed) {
    set.seed(init_seed)
    feats <- mcg_lop_features(Y, order = order)
    km <- tryCatch(stats::kmeans(scale(feats), centers = J, nstart = 1, iter.max = 50),
                   error = function(e) NULL)
    z0 <- if (is.null(km)) sample.int(J, n, replace = TRUE) else km$cluster
    # init params
    p <- as.numeric(table(factor(z0, levels = 1:J)) / n); p[p < 1e-6] <- 1e-6; p <- p / sum(p)
    phi <- if (!is.null(fix_phi)) fix_phi else 0.3
    v_sq <- pmax(apply(Y, 2, var), 1e-4)
    if (v_mode == "single") v_sq <- rep(mean(v_sq), d)
    beta <- matrix(0, J, P)
    for (j in 1:J) {
      idx <- which(z0 == j)
      ybar <- if (length(idx) <= 1) colMeans(Y) else colMeans(Y[idx, , drop = FALSE])
      beta[j, ] <- solve_beta(ybar, phi, v_sq)
    }
    mu <- beta %*% t(basis)

    ll_old <- -Inf
    for (it in seq_len(max_iter)) {
      # E-step: responsibilities
      logcomp <- matrix(0, n, J)
      for (j in 1:J) logcomp[, j] <- log(p[j]) + .mcg_sad_loglik_rows(Y, mu[j, ], phi, v_sq, resets)
      m <- apply(logcomp, 1, max)
      lse <- m + log(rowSums(exp(logcomp - m)))
      ll <- sum(lse)
      gamma <- exp(logcomp - lse)
      if (is.finite(ll) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) { ll_old <- ll; break }
      ll_old <- ll
      # M-step
      Nk <- colSums(gamma); p <- pmax(Nk / n, 1e-8); p <- p / sum(p)
      for (j in 1:J) {
        w <- gamma[, j]; sw <- sum(w)
        ybar <- if (sw < 1e-8) colMeans(Y) else colSums(Y * w) / sw
        beta[j, ] <- solve_beta(ybar, phi, v_sq)
      }
      mu <- beta %*% t(basis)
      # phi (profile) and v
      if (is.null(fix_phi)) {
        negQ <- function(ph) {
          tot <- 0
          ssq <- numeric(d)
          for (j in 1:J) {
            r <- sweep(Y, 2, mu[j, ], "-"); s <- .mcg_whiten(r, ph, resets)
            ssq <- ssq + colSums(gamma[, j] * s^2)
          }
          if (v_mode == "single") {
            vv <- rep(sum(ssq) / (n * d), d)
          } else vv <- pmax(ssq / n, 1e-8)
          0.5 * (sum(n * log(2 * pi * vv)) + sum(ssq / vv))
        }
        opt <- optimize(negQ, c(-0.95, 0.95))
        phi <- opt$minimum
      }
      ssq <- numeric(d)
      for (j in 1:J) {
        r <- sweep(Y, 2, mu[j, ], "-"); s <- .mcg_whiten(r, phi, resets)
        ssq <- ssq + colSums(gamma[, j] * s^2)
      }
      v_sq <- if (v_mode == "single") rep(sum(ssq) / (n * d), d) else pmax(ssq / n, 1e-8)
    }
    z <- max.col(gamma, ties.method = "first")
    list(z = z, loglik = ll_old, p = p, phi = phi, v_sq = v_sq, beta = beta)
  }

  best <- NULL; all_inits <- vector("list", n_init)
  for (k in seq_len(n_init)) {
    fit <- tryCatch(one_run(seed + 7919L * k), error = function(e) NULL)
    all_inits[[k]] <- fit
    if (!is.null(fit) && (is.null(best) || fit$loglik > best$loglik)) best <- fit
    if (verbose) cat(sprintf("  EM init %d loglik=%.1f\n", k, if (is.null(fit)) NA else fit$loglik))
  }
  best$all_inits <- all_inits
  best
}

# =========================================================================
# Generalised SAD(1) MCMC over arbitrary block structure (ablation kernel).
# =========================================================================
# A single nimble distribution handling any block layout via a reset indicator.

.mcg_dSADblocks <- nimble::nimbleFunction(
  run = function(x = double(1), mean = double(1), phi = double(0),
                 v_sq = double(1), reset = double(1), d = integer(0),
                 log = integer(0, default = 0L)) {
    returnType(double(0))
    q <- 0.0; s <- 0.0
    for (t in 1:d) {
      if (reset[t] > 0.5) s <- (x[t] - mean[t]) else s <- (x[t] - mean[t]) + phi * s
      q <- q + (s * s) / v_sq[t] + log(2.0 * 3.141592653589793 * v_sq[t])
    }
    logdens <- -0.5 * q
    if (log) return(logdens) else return(exp(logdens))
  }
)
.mcg_rSADblocks <- nimble::nimbleFunction(
  run = function(n = integer(0, default = 1), mean = double(1), phi = double(0),
                 v_sq = double(1), reset = double(1), d = integer(0)) {
    returnType(double(1))
    out <- numeric(d, init = TRUE); e_prev <- 0.0
    for (t in 1:d) {
      e <- rnorm(1, 0, sqrt(v_sq[t]))
      if (reset[t] > 0.5) e_prev <- 0.0
      out[t] <- mean[t] + e - phi * e_prev
      e_prev <- e
    }
    return(out)
  }
)
.mcg_sadblocks_registered <- new.env(parent = emptyenv())
.mcg_sadblocks_registered$done <- FALSE
mcg_register_sadblocks <- function() {
  assign("dSADblocks", .mcg_dSADblocks, envir = .GlobalEnv)
  assign("rSADblocks", .mcg_rSADblocks, envir = .GlobalEnv)
  if (.mcg_sadblocks_registered$done) return(invisible(TRUE))
  nimble::registerDistributions(list(
    dSADblocks = list(
      BUGSdist = "dSADblocks(mean, phi, v_sq, reset, d)",
      types = c("value=double(1)", "mean=double(1)", "phi=double(0)",
                "v_sq=double(1)", "reset=double(1)", "d=integer(0)"),
      discrete = FALSE)
  ))
  .mcg_sadblocks_registered$done <- TRUE
  invisible(TRUE)
}

# Generalised runner. blocks: list of column-index vectors. basis: d x P design.
# v_mode: "pertime"|"single"; fix_phi: NULL or fixed value (0 = independence).
mcg_run_sad_mcmc <- function(Y, J, blocks, basis, v_mode = c("pertime", "single"),
                             fix_phi = NULL, niter = 8000L, thin = 1L,
                             burnin_frac = 0.25, seed = 1L, init_z = NULL,
                             priors = list(alpha_v = 1, beta_v = 1, mu_phi = 0.25,
                                           eta_phi = 1, sigma_beta = 0.5),
                             trace_draws = 0L) {
  v_mode <- match.arg(v_mode)
  mcg_register_sadblocks()
  set.seed(seed)
  Y <- as.matrix(Y); storage.mode(Y) <- "double"
  n <- nrow(Y); d <- ncol(Y); P <- ncol(basis)
  reset <- numeric(d); for (b in blocks) reset[b[1]] <- 1
  nv <- if (v_mode == "single") 1L else d
  vmap <- if (v_mode == "single") rep(1L, d) else seq_len(d)

  XtX <- crossprod(basis)
  mu_beta <- as.numeric(solve(XtX + diag(1e-6, P), crossprod(basis, colMeans(Y))))
  sigma_beta <- rep(priors$sigma_beta, P)

  # init via kmeans on LOP features, or from a caller-supplied labelling
  # (init_z lets the EM-vs-MCMC init-sensitivity study start both engines from
  # identical diverse starts).
  if (is.null(init_z)) {
    feats <- mcg_lop_features(Y, order = 4L)
    km <- stats::kmeans(scale(feats), centers = J, nstart = 10, iter.max = 100)
    z0 <- km$cluster
  } else {
    z0 <- as.integer(init_z)
  }
  beta0 <- matrix(0, J, P)
  for (j in 1:J) {
    idx <- which(z0 == j); yb <- if (length(idx) <= 1) colMeans(Y) else colMeans(Y[idx, , drop = FALSE])
    beta0[j, ] <- as.numeric(solve(XtX + diag(1e-6, P), crossprod(basis, yb)))
  }
  v0 <- pmax(apply(Y, 2, var), 1e-4)
  p0 <- as.numeric(table(factor(z0, levels = 1:J)) / n)

  # Two model paths. The common case (free phi, per-time variance) is
  # structurally identical to the published run_mcmc_binary model, so it keeps
  # conjugate p sampling and runs at full speed. The variant case (single
  # variance and/or fixed phi, used only by the 1C covariance ablation) uses a
  # general parameterisation with auto-conjugacy disabled.
  common <- is.null(fix_phi) && v_mode == "pertime"

  if (common) {
    code <- nimble::nimbleCode({
      for (j in 1:J) alpha[j] <- 1
      p[1:J] ~ ddirch(alpha[1:J])
      for (j in 1:J) for (k in 1:P) beta[j, k] ~ dnorm(mu_beta[k], tau = 1 / sigma_beta[k]^2)
      for (t in 1:d) v_sq[t] ~ dinvgamma(alpha_v, beta_v)
      phi ~ T(dnorm(mu_phi, tau = 1 / eta_phi), -1, 1)
      for (i in 1:n) z[i] ~ dcat(p[1:J])
      for (j in 1:J) for (t in 1:d) mu[j, t] <- inprod(beta[j, 1:P], Z0[t, 1:P])
      for (i in 1:n) y[i, 1:d] ~ dSADblocks(mean = mu[z[i], 1:d], phi = phi,
                                            v_sq = v_sq[1:d], reset = reset[1:d], d = d)
    })
    constants <- list(n = n, d = d, J = as.integer(J), P = as.integer(P),
                      Z0 = basis, reset = reset,
                      alpha_v = priors$alpha_v, beta_v = priors$beta_v,
                      mu_phi = priors$mu_phi, eta_phi = priors$eta_phi,
                      sigma_beta = sigma_beta, mu_beta = mu_beta)
    inits <- list(z = z0, beta = beta0, v_sq = v0, phi = 0.3, p = p0)
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = list(y = Y), inits = inits, check = FALSE)
    cmodel <- nimble::compileNimble(model, showCompilerOutput = FALSE)
    conf <- nimble::configureMCMC(model, print = FALSE)
    conf$setMonitors(c("phi", "p", "v_sq", "beta")); conf$addMonitors("z")
    conf$removeSamplers("z"); for (i in 1:n) conf$addSampler(paste0("z[", i, "]"), type = "categorical")
    conf$removeSamplers("phi")
    conf$addSampler("phi", type = "slice",
                    control = list(adaptive = TRUE, adaptInterval = 200,
                                   lower = -1 + 1e-8, upper = 1 - 1e-8))
    for (t in 1:d) {
      conf$removeSamplers(paste0("v_sq[", t, "]"))
      conf$addSampler(paste0("v_sq[", t, "]"), type = "slice",
                      control = list(adaptive = TRUE, adaptInterval = 200, lower = 1e-8))
    }
  } else {
    use_phi_node <- is.null(fix_phi); phi_const <- if (use_phi_node) 0 else fix_phi
    code <- nimble::nimbleCode({
      for (j in 1:J) alpha[j] <- 1
      p[1:J] ~ ddirch(alpha[1:J])
      for (j in 1:J) for (k in 1:P) beta[j, k] ~ dnorm(mu_beta[k], tau = 1 / sigma_beta[k]^2)
      for (t in 1:nv) v_raw[t] ~ dinvgamma(alpha_v, beta_v)
      for (t in 1:d) v_sq[t] <- v_raw[vmap[t]]
      phi_free ~ T(dnorm(mu_phi, tau = 1 / eta_phi), -1, 1)
      phi <- use_phi * phi_free + (1 - use_phi) * phi_fixed
      for (i in 1:n) z[i] ~ dcat(p[1:J])
      for (j in 1:J) for (t in 1:d) mu[j, t] <- inprod(beta[j, 1:P], Z0[t, 1:P])
      for (i in 1:n) y[i, 1:d] ~ dSADblocks(mean = mu[z[i], 1:d], phi = phi,
                                            v_sq = v_sq[1:d], reset = reset[1:d], d = d)
    })
    constants <- list(n = n, d = d, J = as.integer(J), P = as.integer(P), nv = nv,
                      vmap = as.integer(vmap), Z0 = basis, reset = reset,
                      alpha_v = priors$alpha_v, beta_v = priors$beta_v,
                      mu_phi = priors$mu_phi, eta_phi = priors$eta_phi,
                      sigma_beta = sigma_beta, mu_beta = mu_beta,
                      use_phi = as.numeric(use_phi_node), phi_fixed = phi_const)
    inits <- list(z = z0, beta = beta0,
                  v_raw = if (v_mode == "single") mean(v0) else v0,
                  phi_free = 0.3, p = p0)
    model <- nimble::nimbleModel(code, constants = constants,
                                 data = list(y = Y), inits = inits, check = FALSE)
    cmodel <- nimble::compileNimble(model, showCompilerOutput = FALSE)
    conf <- nimble::configureMCMC(model, print = FALSE, useConjugacy = FALSE)
    conf$setMonitors(c("phi_free", "p", "v_raw", "beta")); conf$addMonitors("z")
    conf$removeSamplers("z"); for (i in 1:n) conf$addSampler(paste0("z[", i, "]"), type = "categorical")
    if (use_phi_node) {
      conf$removeSamplers("phi_free")
      conf$addSampler("phi_free", type = "slice",
                      control = list(adaptive = TRUE, adaptInterval = 200,
                                     lower = -1 + 1e-8, upper = 1 - 1e-8))
    }
    for (t in 1:nv) {
      conf$removeSamplers(paste0("v_raw[", t, "]"))
      conf$addSampler(paste0("v_raw[", t, "]"), type = "slice",
                      control = list(adaptive = TRUE, adaptInterval = 200, lower = 1e-8))
    }
  }
  for (j in 1:J) {
    for (k in 1:P) conf$removeSamplers(sprintf("beta[%d, %d]", j, k))
    conf$addSampler(sprintf("beta[%d, 1:%d]", j, P), type = "RW_block",
                    control = list(adaptive = TRUE, scale = 0.05, adaptInterval = 200))
  }
  mcmc <- nimble::buildMCMC(conf)
  cmcmc <- nimble::compileNimble(mcmc, project = model, showCompilerOutput = FALSE)
  samples <- nimble::runMCMC(cmcmc, niter = niter, nburnin = 0, thin = thin,
                             progressBar = FALSE, samplesAsCodaMCMC = FALSE)
  coln <- colnames(samples)
  trace_draws <- max(0L, as.integer(trace_draws))
  trace_iterations_full <- integer()
  trace_samples_full <- NULL
  if (trace_draws > 0L) {
    n_trace <- min(trace_draws, nrow(samples))
    trace_rows <- unique(as.integer(round(seq(
      1L, nrow(samples), length.out = n_trace
    ))))
    trace_parameter_columns <- !grepl("^z\\[", coln)
    trace_iterations_full <- trace_rows * thin
    trace_samples_full <- samples[
      trace_rows, trace_parameter_columns, drop = FALSE
    ]
  }
  cut <- floor(nrow(samples) * burnin_frac)
  post <- samples[(cut + 1L):nrow(samples), , drop = FALSE]
  S <- nrow(post)
  z_cols <- grep("^z\\[", coln, value = TRUE)
  z_idx <- as.integer(sub("^z\\[(\\d+)\\]$", "\\1", z_cols))
  z_post <- matrix(as.integer(round(post[, z_cols, drop = FALSE])), nrow = S)[, order(z_idx), drop = FALSE]
  rl <- BPFC:::.relabel_ecr(z_post, J = J)
  z_mode <- apply(rl$z, 2, function(x) which.max(tabulate(x, nbins = J)))
  list(
    clustering = z_mode,
    samples = post,
    z_post = rl$z,
    relabel_perm = rl$perm,
    trace_iterations_full = trace_iterations_full,
    trace_samples_full = trace_samples_full,
    model_info = list(J = J, n = n, d = d)
  )
}

# ---- baseline method wrappers --------------------------------------------

mcg_method_kmeans <- function(Y, K, seed = 1L, order = 4L) {
  set.seed(seed)
  feats <- mcg_lop_features(Y, order = order)
  stats::kmeans(feats, centers = K, nstart = 50, iter.max = 100)$cluster
}

# Supervised oracle baselines: one stratified 70/30 split (seeded) -> metrics
# computed on the held-out test fold. Returns the predicted labels aligned to
# the test indices plus the test truth (caller computes metrics).
mcg_stratified_split <- function(z, p = 0.7, seed = 1L) {
  set.seed(seed)
  tr <- integer(0)
  for (cl in unique(z)) {
    idx <- which(z == cl); k <- max(1L, floor(length(idx) * p))
    tr <- c(tr, sample(idx, k))
  }
  list(train = sort(tr), test = setdiff(seq_along(z), tr))
}

mcg_method_svm <- function(Y, z_true, seed = 1L) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    warning("e1071 not installed; SVM baseline skipped (NA).")
    return(list(pred = NA_integer_, truth = NA_integer_))
  }
  sp <- mcg_stratified_split(z_true, 0.7, seed)
  Ytr <- Y[sp$train, , drop = FALSE]; Yte <- Y[sp$test, , drop = FALSE]
  ztr <- z_true[sp$train]; zte <- z_true[sp$test]
  m <- colMeans(Ytr); s <- apply(Ytr, 2, sd); s[s < 1e-8] <- 1
  Ytr <- scale(Ytr, m, s); Yte <- scale(Yte, m, s)
  fit <- e1071::svm(Ytr, as.factor(ztr), kernel = "radial", scale = FALSE,
                    type = "C-classification")
  pred <- as.integer(as.character(predict(fit, Yte)))
  list(pred = pred, truth = zte)
}

mcg_method_xgb <- function(Y, z_true, seed = 1L, nrounds = 200L) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    warning("xgboost not installed; XGBoost baseline skipped (NA).")
    return(list(pred = NA_integer_, truth = NA_integer_))
  }
  sp <- mcg_stratified_split(z_true, 0.7, seed)
  levs <- sort(unique(z_true))
  Ytr <- Y[sp$train, , drop = FALSE]; Yte <- Y[sp$test, , drop = FALSE]
  ztr <- match(z_true[sp$train], levs) - 1; zte <- z_true[sp$test]
  m <- colMeans(Ytr); s <- apply(Ytr, 2, sd); s[s < 1e-8] <- 1
  Ytr <- scale(Ytr, m, s); Yte <- scale(Yte, m, s)
  dtr <- xgboost::xgb.DMatrix(Ytr, label = ztr)
  fit <- xgboost::xgb.train(list(objective = "multi:softmax", num_class = length(levs),
                                 max_depth = 6, eta = 0.3, subsample = 0.8,
                                 colsample_bytree = 0.8),
                            dtr, nrounds = nrounds, verbose = 0)
  pred <- levs[as.integer(predict(fit, xgboost::xgb.DMatrix(Yte))) + 1]
  list(pred = pred, truth = zte)
}

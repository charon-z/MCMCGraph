#' Run Bayesian MCMC functional clustering for binary (K=2) longitudinal data
#'
#' Input can be wide (n x 2d) or long (id/ind/time/y). K is fixed at 2 and order fixed at 4.
#'
#' @param x data (matrix/data.frame). Wide: n x 2d with concat layout. Long: columns id/ind/time/y.
#' @param J number of clusters
#' @param format "auto" | "wide" | "long"
#' @param times optional numeric vector length d (for basis scaling)
#' @param niter total iterations
#' @param thin thinning
#' @param n_batch batch count to allow checkpoints outside package
#' @param burnin_frac fraction discarded for z_mode (default 0.25)
#' @param priors list of prior settings (optional)
#' @param init optional initialization list with `z`, or complete `z`, `beta`,
#'   `v_sq`, `phi` and `p`. The default `NULL` keeps the original k-means
#'   initialization.
#' @param seed random seed
#' @param two_phi logical. If `FALSE` (default), both states share a single
#'   SAD(1) time-correlation parameter `phi` (the special case phi = psi used
#'   for all reported analyses). If `TRUE`, fit the general model with a
#'   state-specific correlation parameter (`phi` for state 1, `psi` for state 2).
#' @return An object of class `mcmcgraph_result`: a list with `clustering` (MAP
#'   cluster label per feature), `cluster_counts`, `cluster_prob` (n x J
#'   posterior cluster probabilities), `cluster_uncertainty` (1 - max posterior
#'   probability), `posterior_mean` (list of `beta`, `p`, `v_sq`, `phi`, `psi`),
#'   `posterior_samples` (post-burn-in, relabeled draws of `phi`, `psi`, `p`,
#'   `v_sq`, `beta`, `z` and `all_params`; `psi` is `NULL` unless `two_phi`),
#'   `model_info` and `data_info`.
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' fit <- run_mcmc_binary(example_binary$y, J = 3, times = example_binary$times,
#'                        niter = 300, seed = 1)
#' table(fit$clustering)
#' head(fit$cluster_prob)
#' }
run_mcmc_binary <- function(
    x,
    J,
    format = c("auto", "wide", "long"),
    times = NULL,
    niter = 3000,
    thin = 1,
    n_batch = 5L,
    burnin_frac = 0.25,
    priors = list(alpha_v = 1, beta_v = 1, mu_phi = 0.25, eta_phi = 1, sigma_beta = 0.5),
    init = NULL,
    seed = 123,
    two_phi = FALSE
) {
  .mcmcgraph_register_sad()
  format <- match.arg(format)
  if (!requireNamespace("nimble", quietly = TRUE)) stop("Package 'nimble' is required.")

  set.seed(seed)

  # data adapt
  dat <- as_binary_data(x, format = format, times = times)
  y <- dat$y
  d_single <- dat$d_single
  times_single <- dat$times_single
  n <- nrow(y)
  d <- ncol(y)           # 2*d_single

  # basis (order fixed=4 => P=10)
  Z0_binary <- make_Z0_binary(times_single)
  P <- ncol(Z0_binary)   # 10

  # priors
  alpha_v <- priors$alpha_v
  beta_v  <- priors$beta_v
  mu_phi  <- priors$mu_phi
  eta_phi <- priors$eta_phi
  sigma_beta <- rep(priors$sigma_beta, P)

  # mu_beta from global mean curve
  XTX <- crossprod(Z0_binary)
  XTy <- crossprod(Z0_binary, colMeans(y))
  mu_beta <- as.numeric(solve(XTX + diag(1e-6, P), XTy))

  # register distribution (one-time)
  .mcmcgraph_register_sad()

  # model code (P fixed=10, but we keep P for safety).
  # two_phi = FALSE (default): a single time-correlation parameter phi shared by
  #   both states (paper Section 2.3 special case phi = psi); reproduces all
  #   results reported in the manuscript.
  # two_phi = TRUE: the general SAD(1) model with a state-specific correlation
  #   parameter (phi for state 1, psi for state 2) via dSADmvnorm2.
  if (!two_phi) {
    modelCode <- nimble::nimbleCode({
      for (j in 1:J) alpha[j] <- 1
      p[1:J] ~ ddirch(alpha[1:J])

      for (j in 1:J) {
        for (k in 1:P) {
          beta[j, k] ~ dnorm(mu_beta[k], tau = 1 / sigma_beta[k]^2)
        }
      }

      for (t in 1:d) {
        v_sq[t] ~ dinvgamma(alpha_v, beta_v)
      }

      phi ~ T(dnorm(mu_phi, tau = 1 / eta_phi), -1, 1)

      for (i in 1:n) {
        z[i] ~ dcat(p[1:J])
      }

      for (j in 1:J) {
        for (t in 1:d) {
          mu[j, t] <- inprod(beta[j, 1:P], Z0[t, 1:P])
        }
      }

      for (i in 1:n) {
        y[i, 1:d] ~ dSADmvnorm(
          mean = mu[z[i], 1:d],
          phi = phi,
          v_sq = v_sq[1:d],
          d1 = d_single,
          d2 = d_single
        )
      }
    })
  } else {
    modelCode <- nimble::nimbleCode({
      for (j in 1:J) alpha[j] <- 1
      p[1:J] ~ ddirch(alpha[1:J])

      for (j in 1:J) {
        for (k in 1:P) {
          beta[j, k] ~ dnorm(mu_beta[k], tau = 1 / sigma_beta[k]^2)
        }
      }

      for (t in 1:d) {
        v_sq[t] ~ dinvgamma(alpha_v, beta_v)
      }

      # State-specific correlation parameters, both truncated to (-1, 1).
      phi ~ T(dnorm(mu_phi, tau = 1 / eta_phi), -1, 1)
      psi ~ T(dnorm(mu_phi, tau = 1 / eta_phi), -1, 1)

      for (i in 1:n) {
        z[i] ~ dcat(p[1:J])
      }

      for (j in 1:J) {
        for (t in 1:d) {
          mu[j, t] <- inprod(beta[j, 1:P], Z0[t, 1:P])
        }
      }

      for (i in 1:n) {
        y[i, 1:d] ~ dSADmvnorm2(
          mean = mu[z[i], 1:d],
          phi = phi,
          psi = psi,
          v_sq = v_sq[1:d],
          d1 = d_single,
          d2 = d_single
        )
      }
    })
  }

  constants <- list(
    n = n,
    d = d,
    d_single = d_single,
    J = as.integer(J),
    P = as.integer(P),
    Z0 = Z0_binary,
    alpha_v = alpha_v,
    beta_v  = beta_v,
    mu_phi = mu_phi,
    eta_phi = eta_phi,
    sigma_beta = sigma_beta,
    mu_beta = mu_beta
  )

  data_list <- list(y = as.matrix(y))

  # init. NULL preserves the original k-means initialization; callers can pass
  # a neural warm-start or any other data-adaptive initial allocation.
  inits_full <- validate_initialization_binary(init, y, J, Z0_binary)
  init_method <- if (is.null(init)) "kmeans" else if (!is.null(init$method)) init$method else "user"
  inits <- inits_full[c("z", "beta", "v_sq", "phi", "p")]
  if (two_phi) inits$psi <- inits$phi   # initialize state-2 correlation at the state-1 value

  # build model
  model <- nimble::nimbleModel(
    code = modelCode,
    constants = constants,
    data = data_list,
    inits = inits,
    check = FALSE
  )
  cmodel <- nimble::compileNimble(model, showCompilerOutput = FALSE)

  # configure mcmc
  conf <- nimble::configureMCMC(model, print = FALSE)
  conf$setMonitors(c("phi","p","v_sq","beta"))
  conf$addMonitors("z")
  if (two_phi) conf$addMonitors("psi")

  # z samplers
  conf$removeSamplers("z")
  for (i in 1:n) conf$addSampler(target = paste0("z[", i, "]"), type = "categorical")

  # phi slice sampler, bounded to the truncation support (-1, 1)
  conf$removeSamplers("phi")
  conf$addSampler("phi", type = "slice",
                  control = list(adaptive = TRUE, adaptInterval = 200,
                                 lower = -1 + 1e-8, upper = 1 - 1e-8))
  if (two_phi) {
    conf$removeSamplers("psi")
    conf$addSampler("psi", type = "slice",
                    control = list(adaptive = TRUE, adaptInterval = 200,
                                   lower = -1 + 1e-8, upper = 1 - 1e-8))
  }

  # v_sq slice
  for (t in 1:d) {
    conf$removeSamplers(paste0("v_sq[", t, "]"))
    conf$addSampler(paste0("v_sq[", t, "]"),
                    type = "slice",
                    control = list(adaptive = TRUE, adaptInterval = 200, lower = 1e-8))
  }

  # beta block sampler
  for (j in 1:J) {
    for (k in 1:P) conf$removeSamplers(sprintf("beta[%d, %d]", j, k))
    conf$addSampler(
      target = sprintf("beta[%d, 1:%d]", j, P),
      type = "RW_block",
      control = list(adaptive = TRUE, scale = 0.05, adaptInterval = 200)
    )
  }

  mcmc <- nimble::buildMCMC(conf)
  cmcmc <- nimble::compileNimble(mcmc, project = model, showCompilerOutput = FALSE)

  # run in batches (no checkpoint here; caller can do it outside if needed)
  n_batch <- as.integer(max(1L, n_batch))
  iter_per_batch <- ceiling(niter / n_batch)
  samples_list <- vector("list", n_batch)
  curr <- 0L

  # helper: nimble expects list-of-lists for inits if nchains=1 with runMCMC sometimes
  wrap_inits <- function(ini) list(ini)

  inits_next <- inits
  for (b in seq_len(n_batch)) {
    its <- min(iter_per_batch, niter - curr)
    if (its <= 0) break

    sb <- nimble::runMCMC(
      cmcmc,
      niter = its,
      nburnin = 0,
      thin = thin,
      inits = wrap_inits(inits_next),
      progressBar = interactive(),
      samplesAsCodaMCMC = FALSE
    )
    if (inherits(sb, "mcmc")) sb <- as.matrix(sb)
    samples_list[[b]] <- sb
    curr <- curr + its

    # make next inits from last row (lightweight)
    last <- sb[nrow(sb), ]
    inits_next <- .make_inits_from_last(last, J = J, d = d, n = n, P = P)
  }

  samples <- do.call(rbind, samples_list)
  coln <- colnames(samples)

  # ---- discard burn-in (paper 4.1: keep the last (1 - burnin_frac) draws) ----
  burnin_cut <- floor(nrow(samples) * burnin_frac)
  post <- samples[(burnin_cut + 1L):nrow(samples), , drop = FALSE]
  S <- nrow(post)

  # split z vs other parameters
  z_cols <- grep("^z\\[", coln, value = TRUE)
  other_cols <- setdiff(coln, z_cols)

  # canonical, index-ordered arrays (so cluster index j is unambiguous)
  z_post <- matrix(as.integer(round(post[, z_cols, drop = FALSE])), nrow = S)
  # order z columns by index i
  z_idx <- as.integer(sub("^z\\[(\\d+)\\]$", "\\1", z_cols))
  z_post <- z_post[, order(z_idx), drop = FALSE]

  p_cols <- grep("^p\\[", coln, value = TRUE)
  p_idx  <- as.integer(sub("^p\\[(\\d+)\\]$", "\\1", p_cols))
  p_post <- post[, p_cols, drop = FALSE][, order(p_idx), drop = FALSE]  # S x J

  beta_cols <- grep("^beta\\[", coln, value = TRUE)
  bx <- regmatches(beta_cols, regexec("^beta\\[(\\d+),\\s*(\\d+)\\]$", beta_cols))
  bj <- vapply(bx, function(z) as.integer(z[2]), integer(1))
  bk <- vapply(bx, function(z) as.integer(z[3]), integer(1))
  beta_arr <- array(0, dim = c(S, J, P))                # S x J x P
  for (c in seq_along(beta_cols)) beta_arr[, bj[c], bk[c]] <- post[, beta_cols[c]]

  v_cols <- grep("^v_sq\\[", coln, value = TRUE)
  v_idx  <- as.integer(sub("^v_sq\\[(\\d+)\\]$", "\\1", v_cols))
  v_post <- post[, v_cols, drop = FALSE][, order(v_idx), drop = FALSE]  # S x d
  phi_post <- post[, "phi", drop = FALSE]                                # S x 1
  psi_post <- if (two_phi && "psi" %in% coln) post[, "psi", drop = FALSE] else NULL  # S x 1 (general model)

  # ---- resolve label switching, then relabel cluster-indexed draws ----
  rl <- .relabel_ecr(z_post, J = J)
  z_post   <- rl$z
  p_post   <- .apply_perm_vec(p_post, rl$perm)
  beta_arr <- .apply_perm_beta(beta_arr, rl$perm)

  # ---- MAP clustering, posterior probabilities and uncertainty (paper 4.1) ----
  z_mode <- apply(z_post, 2, function(x) which.max(tabulate(x, nbins = J)))
  cluster_prob <- t(apply(z_post, 2, function(x) tabulate(x, nbins = J) / length(x)))
  colnames(cluster_prob) <- paste0("cluster", 1:J)
  z_uncertainty <- 1 - apply(cluster_prob, 1, max)

  # ---- posterior means on the relabeled draws (valid for plug-in BIC) ----
  beta_hat <- apply(beta_arr, c(2, 3), mean)          # J x P
  p_hat    <- colMeans(p_post)
  v_hat    <- colMeans(v_post)
  phi_hat  <- mean(phi_post)
  psi_hat  <- if (!is.null(psi_post)) mean(psi_post) else NA_real_

  res <- list(
    clustering = z_mode,
    cluster_counts = table(z_mode),
    cluster_prob = cluster_prob,
    cluster_uncertainty = z_uncertainty,
    posterior_mean = list(
      beta = beta_hat, p = p_hat, v_sq = v_hat, phi = phi_hat, psi = psi_hat
    ),
    posterior_samples = list(
      phi = phi_post,
      psi = psi_post,
      p = p_post,
      v_sq = v_post,
      beta = beta_arr,
      z = z_post,
      all_params = post[, other_cols, drop = FALSE]
    ),
    model_info = list(
      J = J, n = n, d_single = d_single, d_binary = d,
      order = 4L, P = P,
      niter = niter, thin = thin, burnin_frac = burnin_frac,
      init_method = init_method, two_phi = two_phi
    ),
    data_info = list(
      times_single = times_single,
      col_order = "concat: [i1_t1..i1_td, i2_t1..i2_td]"
    )
  )
  class(res) <- "mcmcgraph_result"
  res
}

#' @keywords internal
.make_inits_from_last <- function(last_row, J, d, n, P) {
  nm <- names(last_row)
  ini <- list()

  if ("phi" %in% nm) ini$phi <- as.numeric(last_row["phi"])
  if ("psi" %in% nm) ini$psi <- as.numeric(last_row["psi"])

  p_names <- grep("^p\\[", nm, value = TRUE)
  if (length(p_names) > 0) {
    p <- numeric(J)
    jidx <- as.integer(sub("^p\\[(\\d+)\\]$", "\\1", p_names))
    p[jidx] <- as.numeric(last_row[p_names])
    p <- pmax(p, 1e-8); p <- p / sum(p)
    ini$p <- p
  }

  v_names <- grep("^v_sq\\[", nm, value = TRUE)
  if (length(v_names) > 0) {
    v <- numeric(d)
    tidx <- as.integer(sub("^v_sq\\[(\\d+)\\]$", "\\1", v_names))
    v[tidx] <- as.numeric(last_row[v_names])
    ini$v_sq <- pmax(v, 1e-8)
  }

  b_names <- grep("^beta\\[", nm, value = TRUE)
  if (length(b_names) > 0) {
    beta <- matrix(0, nrow = J, ncol = P)
    rx <- regexec("^beta\\[(\\d+),\\s*(\\d+)\\]$", b_names)
    parsed <- regmatches(b_names, rx)
    for (i in seq_along(b_names)) {
      jj <- as.integer(parsed[[i]][2]); kk <- as.integer(parsed[[i]][3])
      beta[jj, kk] <- as.numeric(last_row[b_names[i]])
    }
    ini$beta <- beta
  }

  z_names <- grep("^z\\[", nm, value = TRUE)
  if (length(z_names) > 0) {
    z <- integer(n)
    iidx <- as.integer(sub("^z\\[(\\d+)\\]$", "\\1", z_names))
    z[iidx] <- as.integer(round(as.numeric(last_row[z_names])))
    z <- pmin(pmax(z, 1L), as.integer(J))
    ini$z <- z
  }
  ini
}

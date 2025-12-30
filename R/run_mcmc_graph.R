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
#' @param seed random seed
#' @return a list result
#' @export
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
    seed = 123
) {
  format <- match.arg(format)
  if (!requireNamespace("nimble", quietly = TRUE)) stop("Package 'nimble' is required.")
  if (!requireNamespace("coda", quietly = TRUE)) stop("Package 'coda' is required.")

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
  .register_sad_once()


  # model code (P fixed=10, but we keep P for safety)
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

    phi1 ~ dnorm(mu_phi1, tau = 1/eta_phi1)
    phi2 ~ dnorm(mu_phi2, tau = 1/eta_phi2)

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
        phi1 = phi1,
        phi2 = phi2,
        v_sq = v_sq[1:d],
        d1 = d_single,
        d2 = d_single
      )
    }
  })

  constants <- list(
    n = n,
    d = d,
    d_single = d_single,
    J = as.integer(J),
    P = as.integer(P),
    Z0 = Z0_binary,
    alpha_v = alpha_v,
    beta_v  = beta_v,
    mu_phi1 = mu_phi,
    mu_phi2 = mu_phi,
    eta_phi1 = eta_phi,
    eta_phi2 = eta_phi,
    sigma_beta = sigma_beta,
    mu_beta = mu_beta
  )

  data_list <- list(y = as.matrix(y))

  # init
  inits <- initialization_kmeans_binary(y, J, Z0_binary)

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
  conf$setMonitors(c("phi1","phi2","p","v_sq","beta"))
  conf$addMonitors("z")

  # z samplers
  conf$removeSamplers("z")
  for (i in 1:n) conf$addSampler(target = paste0("z[", i, "]"), type = "categorical")

  # phi slice
  conf$removeSamplers("phi1")
  conf$addSampler("phi1", type = "slice", control = list(adaptive = TRUE, adaptInterval = 200))
  conf$removeSamplers("phi2")
  conf$addSampler("phi2", type = "slice", control = list(adaptive = TRUE, adaptInterval = 200))

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

  # split z vs others
  z_cols <- grep("^z\\[", coln, value = TRUE)
  other_cols <- setdiff(coln, z_cols)
  other_samples <- samples[, other_cols, drop = FALSE]

  # z_mode from last (1-burnin_frac) portion
  z_samples <- samples[, z_cols, drop = FALSE]
  burnin_cut <- floor(nrow(z_samples) * burnin_frac)
  z_post <- z_samples[(burnin_cut + 1L):nrow(z_samples), , drop = FALSE]
  z_mode <- apply(z_post, 2, function(x) {
    tb <- table(x)
    as.numeric(names(tb)[which.max(tb)])
  })
  z_uncertainty <- apply(z_post, 2, function(x) {
    tb <- table(x)
    1 - max(tb) / length(x)
  })

  # posterior subsets
  phi_samples <- other_samples[, c("phi1","phi2"), drop = FALSE]
  p_samples   <- other_samples[, grep("^p\\[", colnames(other_samples)), drop = FALSE]
  v_samples   <- other_samples[, grep("^v_sq\\[", colnames(other_samples)), drop = FALSE]
  beta_samples<- other_samples[, grep("^beta\\[", colnames(other_samples)), drop = FALSE]

  res <- list(
    clustering = z_mode,
    cluster_counts = table(z_mode),
    cluster_uncertainty = z_uncertainty,
    posterior_samples = list(
      phi = phi_samples,
      p = p_samples,
      v_sq = v_samples,
      beta = beta_samples,
      all_params = other_samples
    ),
    model_info = list(
      J = J, n = n, d_single = d_single, d_binary = d,
      order = 4L, P = P,
      niter = niter, thin = thin, burnin_frac = burnin_frac
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

  if ("phi1" %in% nm) ini$phi1 <- as.numeric(last_row["phi1"])
  if ("phi2" %in% nm) ini$phi2 <- as.numeric(last_row["phi2"])

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

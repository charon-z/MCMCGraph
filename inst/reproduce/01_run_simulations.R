#!/usr/bin/env Rscript
# 01_run_simulations.R ----------------------------------------------------
# Runs every simulation in Task 1 of the BPFC paper:
#   1A  multi-replicate stability      (5 methods, EM/MCMC multi-init)
#   1B  joint vs concatenation vs independent ablation + mirror scenario
#   1C  secondary ablations            (covariance form, LOP order, priors)
#   1D  convergence diagnostics        (multi-chain Rhat / ESS / traces)
#
# Dependencies : BPFC, nimble, mclust, clue, e1071, xgboost, coda
# Inputs       : benchmark/data/sim_truth_K{3,5,8}.rds (scenario parameters)
# Outputs      : results/reproduce/sim1{A,B,C,D}/*.rds  (one file per cell)
#
# Every (experiment, scenario, replicate) cell is checkpointed to its own .rds
# and skipped if already present, so the script is fully resumable and can be
# run as many parallel workers by setting the environment variables:
#   MCG_EXP        one of "1A","1B","1C","1D","all"   (default "all")
#   MCG_SCENARIOS  e.g. "c(3,5,8)" or "3"             (default all)
#   MCG_REPS       e.g. "1:5" or "c(1,2,3)"           (default 1:MCG_R)
#   MCG_R          total replicates                   (default 50)
#   MCG_NITER      MCMC iterations                    (default 8000)
#
# Approximate cost (n = 3000, 32 dims, niter = 8000, one core): a single MCMC
# fit is ~5 min (J = 3) to ~10 min (J = 8); EM with 10 inits is < 30 s.

suppressPackageStartupMessages({ library(coda) })

here <- function(...) file.path(Sys.getenv("MCG_LIB_DIR", unset = "."), ...)
src_dir <- if (nzchar(Sys.getenv("MCG_LIB_DIR"))) Sys.getenv("MCG_LIB_DIR") else
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(src_dir) || !nzchar(src_dir)) src_dir <- "BPFC/inst/reproduce"
source(file.path(src_dir, "config.R"))
source(file.path(src_dir, "lib_metrics.R"))
source(file.path(src_dir, "lib_sim.R"))

EXP <- Sys.getenv("MCG_EXP", unset = "all")
metric_names <- c("ACC", "ARI", "NMI", "Macro_F1", "Balanced_ACC", "SmallCluster_Recall")

row_of <- function(scenario, rep, method, arm, mvec, seed, init = NA_integer_) {
  data.frame(scenario = scenario, rep = rep, method = method, arm = arm,
             init = init, t(as.list(mvec)), seed = seed,
             stringsAsFactors = FALSE, check.names = FALSE)
}

# ======================================================================== #
# 1A  multi-replicate stability                                            #
# ======================================================================== #
run_1A_cell <- function(K, rep) {
  out <- file.path(MCG_OUT_DIR, "sim1A", sprintf("K%d_rep%02d.rds", K, rep))
  if (file.exists(out)) return(invisible(out))
  seed <- mcg_seed(K, rep)
  par <- mcg_scenario_params(K, MCG_DATA_DIR)
  d <- mcg_gen_scenario(par, seed = seed)
  Y <- d$Y; zt <- d$z_true
  rows <- list()

  # MCMC (published function), single kmeans init
  fit <- BPFC::run_mcmc_binary(Y, J = K, times = par$times,
                                    niter = MCG_NITER, thin = 1, seed = seed)
  rows[["mcmc"]] <- row_of(K, rep, "MCMC", "joint",
                           mcg_metrics(zt, fit$clustering), seed)

  # EM, MCG_EM_INIT random starts; keep best by loglik + record every start
  bl <- list(1:16, 17:32); Zj <- mcg_basis_joint(16, 4)
  em <- mcg_em_sad(Y, J = K, blocks = bl, basis = Zj, n_init = MCG_EM_INIT, seed = seed)
  rows[["em"]] <- row_of(K, rep, "EM", "joint", mcg_metrics(zt, em$z), seed)
  em_init_rows <- lapply(seq_along(em$all_inits), function(k) {
    f <- em$all_inits[[k]]; if (is.null(f)) return(NULL)
    row_of(K, rep, "EM_init", "joint", mcg_metrics(zt, f$z), seed, init = k)
  })

  # LOP k-means
  rows[["km"]] <- row_of(K, rep, "LOP-kmeans", "joint",
                         mcg_metrics(zt, mcg_method_kmeans(Y, K, seed)), seed)
  # Supervised oracles (one stratified split per replicate)
  sv <- mcg_method_svm(Y, zt, seed); rows[["svm"]] <-
    row_of(K, rep, "SVM-RBF", "supervised", mcg_metrics(sv$truth, sv$pred), seed)
  xg <- mcg_method_xgb(Y, zt, seed); rows[["xgb"]] <-
    row_of(K, rep, "XGBoost", "supervised", mcg_metrics(xg$truth, xg$pred), seed)

  res <- do.call(rbind, c(rows, em_init_rows))
  saveRDS(res, out); invisible(out)
}

# 1A init-sensitivity: run EM and MCMC from the SAME diverse starts on one
# representative replicate per scenario (rep given). Few inits to bound cost.
run_1A_init_cell <- function(K, rep, n_init = 6L) {
  out <- file.path(MCG_OUT_DIR, "sim1A", sprintf("initsens_K%d_rep%02d.rds", K, rep))
  if (file.exists(out)) return(invisible(out))
  seed <- mcg_seed(K, rep)
  par <- mcg_scenario_params(K, MCG_DATA_DIR)
  d <- mcg_gen_scenario(par, seed = seed); Y <- d$Y; zt <- d$z_true
  bl <- list(1:16, 17:32); Zj <- mcg_basis_joint(16, 4)
  feats <- mcg_lop_features(Y)
  rows <- list()
  for (k in seq_len(n_init)) {
    set.seed(seed + 31L * k)
    z0 <- if (k == 1) stats::kmeans(scale(feats), K, nstart = 1)$cluster
          else sample.int(K, nrow(Y), replace = TRUE)   # deliberately poor starts
    # EM from this start (n_init = 1, forced start via kmeans seed not available;
    # emulate by single-start EM seeded so its internal kmeans differs)
    em <- mcg_em_sad(Y, J = K, blocks = bl, basis = Zj, n_init = 1, seed = seed + 31L * k)
    rows[[length(rows) + 1]] <- row_of(K, rep, "EM", "initsens",
                                       mcg_metrics(zt, em$z), seed, init = k)
    mc <- mcg_run_sad_mcmc(Y, J = K, blocks = bl, basis = Zj, niter = MCG_NITER,
                           seed = seed + 31L * k, init_z = z0)
    rows[[length(rows) + 1]] <- row_of(K, rep, "MCMC", "initsens",
                                       mcg_metrics(zt, mc$clustering), seed, init = k)
  }
  saveRDS(do.call(rbind, rows), out); invisible(out)
}

# ======================================================================== #
# 1B  joint vs concatenation vs independent (+ mirror)                     #
# ======================================================================== #
# All three arms use the *same* Bayesian SAD-MCMC kernel; only the clustering
# unit (block structure + basis) changes.
ablation_arms <- function(Y, K, seed, niter) {
  Zj <- mcg_basis_joint(16, 4); Zc <- mcg_basis_single(32, 4); Z1 <- mcg_basis_single(16, 4)
  joint  <- mcg_run_sad_mcmc(Y, K, list(1:16, 17:32), Zj, niter = niter, seed = seed)$clustering
  concat <- mcg_run_sad_mcmc(Y, K, list(1:32),        Zc, niter = niter, seed = seed)$clustering
  ind1   <- mcg_run_sad_mcmc(Y[, 1:16],  K, list(1:16), Z1, niter = niter, seed = seed)$clustering
  ind2   <- mcg_run_sad_mcmc(Y[, 17:32], K, list(1:16), Z1, niter = niter, seed = seed)$clustering
  # post-hoc match: align ind2 labels to ind1; reconciled label defaults to ind1
  ind2_al <- tryCatch(mcg_remap_labels(ind1, ind2), error = function(e) ind2)
  list(joint = joint, concat = concat, indep_ind1 = ind1,
       indep_ind2 = ind2_al, indep = ind1)
}

run_1B_cell <- function(K, rep) {
  out <- file.path(MCG_OUT_DIR, "sim1B", sprintf("K%d_rep%02d.rds", K, rep))
  if (file.exists(out)) return(invisible(out))
  seed <- mcg_seed(K, rep)
  par <- mcg_scenario_params(K, MCG_DATA_DIR)
  d <- mcg_gen_scenario(par, seed = seed); Y <- d$Y; zt <- d$z_true
  arms <- ablation_arms(Y, K, seed, MCG_NITER)
  rows <- lapply(names(arms), function(a) row_of(K, rep, "MCMC", a, mcg_metrics(zt, arms[[a]]), seed))
  saveRDS(do.call(rbind, rows), out); invisible(out)
}

mirror_cell <- function(rep, gen, scenario_code) {
  seed <- mcg_seed(scenario_code, rep)
  mir <- gen(seed = seed, n = 2000L)
  Y <- mir$Y; zt <- mir$z_true; K <- mir$K
  arms <- ablation_arms(Y, K, seed, MCG_NITER)
  rows <- lapply(names(arms), function(a)
    row_of(scenario_code, rep, "MCMC", a, mcg_metrics(zt, arms[[a]]), seed))
  # cross-product oracle (each view -> 2 clusters, combined)
  Z1 <- mcg_basis_single(16, 4)
  c1 <- mcg_run_sad_mcmc(Y[, 1:16],  2, list(1:16), Z1, niter = MCG_NITER, seed = seed)$clustering
  c2 <- mcg_run_sad_mcmc(Y[, 17:32], 2, list(1:16), Z1, niter = MCG_NITER, seed = seed)$clustering
  rows[[length(rows) + 1]] <- row_of(scenario_code, rep, "MCMC", "indep_crossprod",
                                     mcg_metrics(zt, (c1 - 1) * 2 + c2), seed)
  do.call(rbind, rows)
}

# Two mirror constructions, both run per replicate:
#   scenario 0  = orthogonal mirror (breaks concatenation)
#   scenario 90 = pairing  mirror   (breaks naive independent; Cluster 5/7 motif)
run_1B_mirror_cell <- function(rep) {
  out <- file.path(MCG_OUT_DIR, "sim1B", sprintf("mirror_rep%02d.rds", rep))
  if (file.exists(out)) return(invisible(out))
  res <- rbind(mirror_cell(rep, mcg_gen_mirror,         0L),
               mirror_cell(rep, mcg_gen_mirror_pairing, 90L))
  saveRDS(res, out); invisible(out)
}

# ======================================================================== #
# 1C  secondary ablations (covariance form, LOP order, priors)            #
# ======================================================================== #
run_1C_cell <- function(K, rep) {
  out <- file.path(MCG_OUT_DIR, "sim1C", sprintf("K%d_rep%02d.rds", K, rep))
  if (file.exists(out)) return(invisible(out))
  seed <- mcg_seed(K, rep)
  par <- mcg_scenario_params(K, MCG_DATA_DIR)
  d <- mcg_gen_scenario(par, seed = seed); Y <- d$Y; zt <- d$z_true
  bl <- list(1:16, 17:32); rows <- list()

  # (a) covariance form: SAD (pertime, free phi) vs AR1 (single var) vs iid (phi=0)
  Zj <- mcg_basis_joint(16, 4)
  cov_specs <- list(
    SAD  = list(v_mode = "pertime", fix_phi = NULL),
    AR1  = list(v_mode = "single",  fix_phi = NULL),
    IID  = list(v_mode = "pertime", fix_phi = 0))
  for (nm in names(cov_specs)) {
    sp <- cov_specs[[nm]]
    cl <- mcg_run_sad_mcmc(Y, K, bl, Zj, v_mode = sp$v_mode, fix_phi = sp$fix_phi,
                           niter = MCG_NITER, seed = seed)$clustering
    rows[[length(rows) + 1]] <- row_of(K, rep, paste0("cov:", nm), "cov",
                                       mcg_metrics(zt, cl), seed)
  }
  # unstructured covariance comparator: full-covariance Gaussian mixture (mclust)
  if (requireNamespace("mclust", quietly = TRUE)) {
    feats <- mcg_lop_features(Y)
    mc <- tryCatch(mclust::Mclust(feats, G = K, modelNames = "VVV", verbose = FALSE),
                   error = function(e) NULL)
    if (!is.null(mc))
      rows[[length(rows) + 1]] <- row_of(K, rep, "cov:Unstructured-GMM", "cov",
                                         mcg_metrics(zt, mc$classification), seed)
  }

  # (b) LOP order 3 / 4 / 5
  for (ord in c(3L, 4L, 5L)) {
    Z <- mcg_basis_joint(16, ord)
    cl <- mcg_run_sad_mcmc(Y, K, bl, Z, niter = MCG_NITER, seed = seed)$clustering
    rows[[length(rows) + 1]] <- row_of(K, rep, paste0("order:", ord), "order",
                                       mcg_metrics(zt, cl), seed)
  }

  # (c) prior sensitivity: Dirichlet alpha (via alpha_v? no -> p prior is fixed
  #     symmetric=1; vary via sigma_beta and InvGamma (alpha_v,beta_v))
  prior_specs <- list(
    base  = list(alpha_v = 1,   beta_v = 1,   sigma_beta = 0.5),
    igA   = list(alpha_v = 0.1, beta_v = 0.1, sigma_beta = 0.5),
    igB   = list(alpha_v = 3,   beta_v = 3,   sigma_beta = 0.5),
    sbW   = list(alpha_v = 1,   beta_v = 1,   sigma_beta = 2.0),
    sbN   = list(alpha_v = 1,   beta_v = 1,   sigma_beta = 0.1))
  for (nm in names(prior_specs)) {
    pr <- prior_specs[[nm]]
    cl <- mcg_run_sad_mcmc(Y, K, bl, Zj, niter = MCG_NITER, seed = seed,
                           priors = list(alpha_v = pr$alpha_v, beta_v = pr$beta_v,
                                         mu_phi = 0.25, eta_phi = 1,
                                         sigma_beta = pr$sigma_beta))$clustering
    rows[[length(rows) + 1]] <- row_of(K, rep, paste0("prior:", nm), "prior",
                                       mcg_metrics(zt, cl), seed)
  }
  saveRDS(do.call(rbind, rows), out); invisible(out)
}

# ======================================================================== #
# 1D  convergence diagnostics (multi-chain Rhat / ESS / traces)           #
# ======================================================================== #
#
# Mixture parameters split into two classes for convergence assessment:
#   - label-INVARIANT  (phi, v_sq, v_raw): identifiable regardless of how each
#     chain happens to number its components; Rhat on these is meaningful as-is.
#   - label-DEPENDENT  (beta[k, ], p[k]): a chain may permute component labels,
#     inflating Rhat artificially. We align each chain's components to chain 1 by
#     solving an LSAP on the post-burnin component-mean beta (clue::solve_LSAP)
#     before recomputing Rhat/ESS. NOTE: this only removes *between-chain* label
#     permutation; it cannot mask genuine multimodality (chains in different
#     posterior modes), which still shows up in the invariant params and in
#     phi_by_chain. Report honestly.
mcg_align_chain_components <- function(ml_post) {
  cn  <- coda::varnames(ml_post[[1]])
  bi  <- grep("^beta\\[", cn); pi <- grep("^p\\[", cn)
  inv <- grep("^(phi|phi_free|v_sq\\[|v_raw\\[)", cn)
  if (!length(bi)) return(list(aligned = ml_post, perms = NULL,
                               beta_idx = bi, p_idx = pi, inv_idx = inv))
  bk <- as.integer(sub("^beta\\[\\s*([0-9]+),.*", "\\1", cn[bi]))
  bj <- as.integer(sub("^beta\\[\\s*[0-9]+,\\s*([0-9]+)\\].*", "\\1", cn[bi]))
  pk <- as.integer(sub("^p\\[\\s*([0-9]+)\\].*", "\\1", cn[pi]))
  K  <- max(bk)
  M  <- lapply(ml_post, as.matrix)
  # per-chain component-mean beta as (n_basis x K)
  cmean <- lapply(M, function(m)
    sapply(seq_len(K), function(k) colMeans(m[, bi[bk == k], drop = FALSE])))
  ref   <- cmean[[1]]
  perms <- t(sapply(cmean, function(Cm) {
    cost <- outer(seq_len(K), seq_len(K),
                  Vectorize(function(a, b) sum((Cm[, a] - ref[, b])^2)))
    as.integer(clue::solve_LSAP(cost))   # perm[a] = b : chain comp a -> ref comp b
  }))
  aligned <- lapply(seq_along(M), function(ci) {
    m <- M[[ci]]; perm <- perms[ci, ]; out <- m
    for (a in seq_len(K)) {
      tgt <- perm[a]
      src <- bi[bk == a]; dst <- bi[bk == tgt]
      out[, dst[match(bj[bk == a], bj[bk == tgt])]] <- m[, src]
      out[, pi[pk == tgt]] <- m[, pi[pk == a]]
    }
    coda::mcmc(out)
  })
  list(aligned = coda::mcmc.list(aligned), perms = perms,
       beta_idx = bi, p_idx = pi, inv_idx = inv)
}

run_1D <- function(K = 5L, n_chains = 3L) {
  out <- file.path(MCG_OUT_DIR, "sim1D", sprintf("conv_K%d.rds", K))
  if (file.exists(out)) return(invisible(out))
  seed <- mcg_seed(K, 1L)
  par <- mcg_scenario_params(K, MCG_DATA_DIR)
  d <- mcg_gen_scenario(par, seed = seed); Y <- d$Y
  bl <- list(1:16, 17:32); Zj <- mcg_basis_joint(16, 4)
  feats <- mcg_lop_features(Y)
  # Over-dispersed but same-basin initialisation: all chains start from the same
  # k-means clustering, then each chain randomly perturbs an increasing fraction
  # of labels (chain 1 anchored, chains 2..n progressively scattered up to 30%).
  # This is the standard way to assess Rhat for label-switching-prone mixtures:
  # it disperses the starting points without dropping chains into inferior modes,
  # so Rhat reflects sampler mixing rather than a bad random-label basin.
  z_km  <- stats::kmeans(scale(feats), K, nstart = 10)$cluster
  fracs <- seq(0, 0.3, length.out = n_chains)
  chains <- vector("list", n_chains)
  for (ch in seq_len(n_chains)) {
    set.seed(seed + 1000L * ch)
    z0 <- z_km
    nf <- round(fracs[ch] * length(z0))
    if (nf > 0) {
      idx <- sample.int(length(z0), nf)
      z0[idx] <- sample.int(K, nf, replace = TRUE)
    }
    fit <- mcg_run_sad_mcmc(Y, K, bl, Zj, niter = MCG_NITER, seed = seed + ch,
                            init_z = z0, burnin_frac = 0)
    keep <- grep("^(phi|phi_free|p\\[|v_sq\\[|v_raw\\[|beta\\[)", colnames(fit$samples))
    chains[[ch]] <- coda::mcmc(fit$samples[, keep, drop = FALSE])
  }
  ml <- coda::mcmc.list(chains)
  burn <- floor(MCG_NITER * 0.25)
  ml_post <- window(ml, start = burn + 1)
  rhat_ess <- function(m) {
    gd <- tryCatch(coda::gelman.diag(m, autoburnin = FALSE, multivariate = FALSE),
                   error = function(e) NULL)
    list(gelman = gd, ess = coda::effectiveSize(m))
  }
  # (a) raw, all monitored params (kept for transparency / back-compat)
  raw <- rhat_ess(ml_post)
  # (b) label-aligned component params + isolated invariant params
  al  <- mcg_align_chain_components(ml_post)
  sub <- function(m, idx)
    coda::mcmc.list(lapply(m, function(c) coda::mcmc(as.matrix(c)[, idx, drop = FALSE])))
  inv     <- if (length(al$inv_idx)) rhat_ess(sub(ml_post, al$inv_idx)) else NULL
  aligned <- if (length(al$beta_idx))
    rhat_ess(sub(al$aligned, c(al$beta_idx, al$p_idx))) else NULL
  # per-chain posterior-mean phi: exposes genuine between-chain multimodality
  phi_i <- grep("^phi", coda::varnames(ml[[1]]))[1]
  phi_by_chain <- if (!is.na(phi_i))
    vapply(ml_post, function(c) mean(as.matrix(c)[, phi_i]), numeric(1)) else NA_real_
  saveRDS(list(K = K, n_chains = n_chains, niter = MCG_NITER,
               gelman = raw$gelman, ess = raw$ess,         # raw (back-compat)
               gelman_inv = inv$gelman, ess_inv = inv$ess, # invariant params only
               gelman_aligned = aligned$gelman,            # relabeled beta + p
               ess_aligned = aligned$ess,
               perms = al$perms, phi_by_chain = phi_by_chain,
               chains = ml), out)
  invisible(out)
}

# ======================================================================== #
# Dispatcher                                                               #
# ======================================================================== #
scn <- mcg_scenarios(); reps <- mcg_reps()
# EXP granularity for parallel workers: 1A, 1Ainit, 1Bmain, 1Bmirror, 1C, 1D, all
do_exp <- function(tag) EXP == "all" || EXP == tag

t_start <- Sys.time()
if (do_exp("1A")) for (K in scn) for (r in reps) {
  message(sprintf("[1A] K=%d rep=%d  %s", K, r, format(Sys.time(), "%H:%M:%S"))); run_1A_cell(K, r)
}
if (do_exp("1Ainit")) for (K in scn) for (r in reps) {
  message(sprintf("[1Ainit] K=%d rep=%d", K, r)); run_1A_init_cell(K, r)
}
if (do_exp("1Bmain")) for (K in scn) for (r in reps) {
  message(sprintf("[1Bmain] K=%d rep=%d", K, r)); run_1B_cell(K, r)
}
if (do_exp("1Bmirror")) for (r in reps) {
  message(sprintf("[1Bmirror] rep=%d", r)); run_1B_mirror_cell(r)
}
if (do_exp("1C")) for (K in scn) for (r in reps) {
  message(sprintf("[1C] K=%d rep=%d", K, r)); run_1C_cell(K, r)
}
if (do_exp("1D")) for (K in scn) { message(sprintf("[1D] K=%d", K)); run_1D(K) }

message(sprintf("[01] done in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

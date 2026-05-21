# relabel.R
# Post-processing relabeling to resolve label switching (paper Section 4.2).
#
# Finite mixture likelihoods are invariant to permutations of the cluster
# labels, so the raw MCMC draws of z, p and beta may switch cluster identities
# between iterations. Averaging such draws directly is meaningless. Here we run
# an ECR-style (Equivalence Classes Representatives) alignment: every iteration
# is relabeled by the permutation that best matches a common reference (pivot)
# allocation, and the pivot is refined until stable.

#' @keywords internal
.all_perms <- function(v) {
  if (length(v) <= 1L) return(list(v))
  out <- list()
  for (i in seq_along(v)) {
    rest <- .all_perms(v[-i])
    for (r in rest) out[[length(out) + 1L]] <- c(v[i], r)
  }
  out
}

#' @keywords internal
#' Best permutation tau (old label j -> tau[j]) maximizing label agreement.
#' Uses the contingency matrix C[j, l] = #{i : z[i] == j and pivot[i] == l};
#' agreement(tau) = sum_j C[j, tau[j]]. When the full permutation matrix PM
#' (n_perm x J) is supplied the optimum is found by a vectorized exhaustive
#' search; otherwise a greedy assignment is used (large J).
.best_perm <- function(z_row, pivot, J, PM = NULL) {
  C <- table(factor(z_row, levels = 1:J), factor(pivot, levels = 1:J))
  C <- matrix(as.numeric(C), nrow = J, ncol = J)

  if (!is.null(PM)) {
    scores <- numeric(nrow(PM))
    for (j in 1:J) scores <- scores + C[j, PM[, j]]   # sum_j C[j, tau_j], vectorized
    return(PM[which.max(scores), ])
  }

  # greedy fallback (J too large to enumerate): assign highest counts first
  tau <- integer(J); used <- logical(J)
  ord <- order(as.numeric(C), decreasing = TRUE)
  for (idx in ord) {
    j <- ((idx - 1L) %% J) + 1L
    l <- ((idx - 1L) %/% J) + 1L
    if (tau[j] == 0L && !used[l]) { tau[j] <- l; used[l] <- TRUE }
  }
  miss_j <- which(tau == 0L); miss_l <- which(!used)
  if (length(miss_j)) tau[miss_j] <- miss_l
  tau
}

#' @keywords internal
#' ECR relabeling. Returns the per-iteration permutations and the relabeled z.
#' @param z_post S x n integer matrix of post-burnin cluster labels.
#' @param J number of clusters.
#' @param max_iter pivot-refinement iterations.
.relabel_ecr <- function(z_post, J, max_iter = 15L) {
  S <- nrow(z_post); n <- ncol(z_post)
  PM <- if (J <= 8L) do.call(rbind, .all_perms(1:J)) else NULL  # enumerate when feasible

  col_mode <- function(M) {
    apply(M, 2, function(col) {
      tb <- tabulate(col, nbins = J)
      which.max(tb)
    })
  }

  pivot <- col_mode(z_post)
  perm_mat <- matrix(rep(1:J, each = S), nrow = S)  # identity by default

  for (it in seq_len(max_iter)) {
    relabeled <- z_post
    for (s in seq_len(S)) {
      tau <- .best_perm(z_post[s, ], pivot, J, PM)
      perm_mat[s, ] <- tau
      relabeled[s, ] <- tau[z_post[s, ]]
    }
    new_pivot <- col_mode(relabeled)
    if (identical(as.integer(new_pivot), as.integer(pivot))) break
    pivot <- new_pivot
  }

  list(perm = perm_mat, z = relabeled, pivot = pivot)
}

#' @keywords internal
#' Apply per-iteration permutations to cluster-indexed draws.
#' tau maps old label j -> new label tau[j]; new[tau[j]] = old[j].
.apply_perm_vec <- function(p_post, perm_mat) {
  S <- nrow(p_post); J <- ncol(p_post)
  out <- p_post
  for (s in seq_len(S)) out[s, perm_mat[s, ]] <- p_post[s, ]
  out
}

#' @keywords internal
#' beta_arr: S x J x P array. Relabel the J dimension per iteration.
.apply_perm_beta <- function(beta_arr, perm_mat) {
  S <- dim(beta_arr)[1]
  out <- beta_arr
  for (s in seq_len(S)) out[s, perm_mat[s, ], ] <- beta_arr[s, , ]
  out
}

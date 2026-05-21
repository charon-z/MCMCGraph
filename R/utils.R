#' @keywords internal
.check_no_na <- function(y) {
  if (anyNA(y)) stop("Input y contains NA. Please impute or remove missing values before running.")
  invisible(TRUE)
}

#' Coerce user data into required binary wide format: n x (2d)
#'
#' Supported:
#' - wide matrix/data.frame: n x (2d), ordered as ind1(t1..td), ind2(t1..td)
#' - long data.frame: columns for id, ind, time, value (ind must have exactly 2 levels)
#'
#' @param x matrix/data.frame
#' @param format "auto" | "wide" | "long"
#' @param times optional numeric vector length d
#' @param id_col, ind_col, time_col, value_col used when long format
#' @return list(y, d_single, times_single)
#' @export
as_binary_data <- function(
    x,
    format = c("auto", "wide", "long"),
    times = NULL,
    id_col = "snp", ind_col = "ind", time_col = "time", value_col = "y"
) {
  format <- match.arg(format)

  if (format == "auto") {
    if (is.data.frame(x) && all(c(ind_col, time_col, value_col) %in% names(x))) {
      format <- "long"
    } else {
      format <- "wide"
    }
  }

  if (format == "wide") {
    y <- as.matrix(x)
    if (!is.numeric(y)) storage.mode(y) <- "double"
    if (ncol(y) %% 2L != 0L) stop("Wide input must have even number of columns = 2*d.")
    d_single <- ncol(y) / 2L
    if (is.null(times)) times <- seq_len(d_single)
    if (length(times) != d_single) stop("times length must equal d_single.")
    .check_no_na(y)
    return(list(y = y, d_single = as.integer(d_single), times_single = as.numeric(times)))
  }

  # long format
  df <- as.data.frame(x)
  need <- c(ind_col, time_col, value_col)
  if (!all(need %in% names(df))) stop("Long input must contain columns: ", paste(need, collapse = ", "))
  inds <- unique(df[[ind_col]])
  if (length(inds) != 2L) stop("Long input must have exactly 2 individuals in column ", ind_col)

  if (is.null(times)) times <- sort(unique(df[[time_col]]))
  times <- as.numeric(times)
  d_single <- length(times)

  # stable ind order
  ind_levels <- sort(as.character(inds))
  df[[ind_col]] <- factor(df[[ind_col]], levels = ind_levels)

  if (!id_col %in% names(df)) stop("Long input must contain an id column: ", id_col)

  ids <- unique(df[[id_col]])
  n <- length(ids)
  y <- matrix(NA_real_, nrow = n, ncol = 2L * d_single)
  rownames(y) <- as.character(ids)

  for (ii in seq_along(ids)) {
    sub <- df[df[[id_col]] == ids[ii], , drop = FALSE]
    for (k in 1:2) {
      ind_k <- levels(df[[ind_col]])[k]
      subk <- sub[sub[[ind_col]] == ind_k, , drop = FALSE]
      vals <- rep(NA_real_, d_single)
      names(vals) <- as.character(times)
      tt <- as.character(subk[[time_col]])
      vals[tt] <- as.numeric(subk[[value_col]])
      y[ii, ((k - 1L) * d_single + 1L):(k * d_single)] <- vals
    }
  }

  .check_no_na(y)
  list(y = y, d_single = as.integer(d_single), times_single = times)
}

# ---- Basis functions ----

#' @keywords internal
calc_Z_legendre <- function(times, order = 4L) {
  if (order != 4L) stop("This package fixes order=4 for the paper setting.")
  d <- length(times)
  q <- 5L
  Z <- matrix(0, d, q)

  tmin <- min(times); tmax <- max(times)
  tx <- if (tmax == tmin) rep(0, d) else (-1 + 2 * (times - tmin) / (tmax - tmin))

  for (i in seq_len(d)) {
    x <- tx[i]
    Z[i, ] <- c(
      1.0,
      x,
      0.5 * (3 * x^2 - 1),
      0.5 * (5 * x^3 - 3 * x),
      0.125 * (35 * x^4 - 30 * x^2 + 3)
    )
  }
  Z
}

#' @keywords internal
make_Z0_binary <- function(times_single) {
  Z1 <- calc_Z_legendre(times_single, order = 4L)  # d x 5
  d <- nrow(Z1); q <- ncol(Z1) # q=5
  Zb <- matrix(0, nrow = 2L * d, ncol = 2L * q)   # (2d) x 10
  Zb[1:d, 1:q] <- Z1
  Zb[(d + 1L):(2L * d), (q + 1L):(2L * q)] <- Z1
  Zb
}

# ---- phi estimation for init ----

#' @keywords internal
estimate_phi_block <- function(mat) {
  n <- nrow(mat); d <- ncol(mat)
  if (d < 2L) return(0)
  yc <- scale(mat, center = TRUE, scale = FALSE)
  num <- 0; den <- 0
  for (i in 1:n) {
    for (t in 1:(d - 1L)) {
      num <- num + yc[i, t] * yc[i, t + 1L]
      den <- den + yc[i, t]^2
    }
  }
  phi <- if (den > 0) num / den else 0
  phi <- max(min(phi, 0.99), -0.99)
  phi
}

# ---- Kmeans init (binary) ----

#' @keywords internal
initialization_kmeans_binary <- function(y, J, Z0_binary) {
  set.seed(123)
  y_scaled <- scale(y)
  km <- stats::kmeans(y_scaled, centers = J, nstart = 30, iter.max = 100)
  z <- km$cluster

  P <- ncol(Z0_binary)  # 10
  beta_init <- matrix(0, nrow = J, ncol = P)
  v_sq_init <- rep(1.0, ncol(y))

  XTX <- crossprod(Z0_binary)  # P x P

  for (j in 1:J) {
    idx <- which(z == j)
    if (length(idx) == 0) next
    sub <- y[idx, , drop = FALSE]
    mean_curve <- if (nrow(sub) == 1) as.numeric(sub) else colMeans(sub)
    XTy <- crossprod(Z0_binary, mean_curve)  # P
    beta_init[j, ] <- as.numeric(solve(XTX + diag(1e-6, P), XTy))
  }

  # residual var init
  # (simple: per timepoint var from residuals of assigned cluster mean)
  res <- matrix(0, nrow = 0, ncol = ncol(y))
  for (j in 1:J) {
    idx <- which(z == j)
    if (length(idx) == 0) next
    fitted <- as.numeric(Z0_binary %*% beta_init[j, ])
    res <- rbind(res, y[idx, , drop = FALSE] - matrix(fitted, nrow = length(idx), ncol = ncol(y), byrow = TRUE))
  }
  if (nrow(res) > 1) {
    v_sq_init <- pmax(apply(res, 2, stats::var), 1e-6)
  }

  d_single <- ncol(y) / 2L
  phi1 <- estimate_phi_block(y[, 1:d_single, drop = FALSE])
  phi2 <- estimate_phi_block(y[, (d_single + 1L):(2L * d_single), drop = FALSE])
  # both individuals share a single phi (paper 2.3): start from the block average,
  # clamped strictly inside the (-1, 1) truncation support
  phi <- max(min((phi1 + phi2) / 2, 0.95), -0.95)

  p_init <- as.numeric(table(z) / length(z))
  list(z = z, beta = beta_init, v_sq = v_sq_init, phi = phi, p = p_init)
}

# ---- log-sum-exp ----
#' @keywords internal
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

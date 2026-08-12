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
#' @param id_col name of the feature-id column (long format only).
#' @param ind_col name of the individual column; must have exactly two levels
#'   (long format only).
#' @param time_col name of the time column (long format only).
#' @param value_col name of the value column (long format only).
#' @return A list with `y` (the n x 2d wide numeric matrix), `d_single` (time
#'   points per individual) and `times_single` (the numeric time grid).
#' @export
#' @examples
#' data(example_binary)
#' dat <- as_binary_data(example_binary$y, format = "wide",
#'                       times = example_binary$times)
#' dim(dat$y)
#' dat$d_single
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
initialization_from_labels_binary <- function(y, z, J, Z0_binary) {
  z <- as.integer(z)
  if (length(z) != nrow(y)) stop("init labels must have length nrow(y).", call. = FALSE)
  if (any(!is.finite(z)) || any(z < 1L) || any(z > J)) {
    stop("init labels must be finite integers in 1:J.", call. = FALSE)
  }
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
  if (length(p_init) < J) p_init <- as.numeric(table(factor(z, levels = 1:J)) / length(z))
  p_init <- pmax(p_init, 1e-8)
  p_init <- p_init / sum(p_init)
  list(z = z, beta = beta_init, v_sq = v_sq_init, phi = phi, p = p_init)
}

#' @keywords internal
initialization_kmeans_binary <- function(y, J, Z0_binary) {
  set.seed(123)
  y_scaled <- scale(y)
  km <- stats::kmeans(y_scaled, centers = J, nstart = 30, iter.max = 100)
  initialization_from_labels_binary(y, km$cluster, J, Z0_binary)
}

#' @keywords internal
validate_initialization_binary <- function(init, y, J, Z0_binary) {
  if (is.null(init)) return(initialization_kmeans_binary(y, J, Z0_binary))
  if (!is.list(init)) stop("init must be a list or NULL.", call. = FALSE)
  if (!is.null(init$z) && (is.null(init$beta) || is.null(init$v_sq) || is.null(init$phi) || is.null(init$p))) {
    init <- initialization_from_labels_binary(y, init$z, J, Z0_binary)
  }

  n <- nrow(y); d <- ncol(y); P <- ncol(Z0_binary)
  need <- c("z", "beta", "v_sq", "phi", "p")
  missing <- setdiff(need, names(init))
  if (length(missing)) stop("init is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  init$z <- as.integer(init$z)
  if (length(init$z) != n || any(init$z < 1L) || any(init$z > J)) {
    stop("init$z must have length nrow(y) and values in 1:J.", call. = FALSE)
  }
  init$beta <- as.matrix(init$beta)
  if (!identical(dim(init$beta), c(as.integer(J), as.integer(P)))) {
    stop("init$beta must be a J x P matrix.", call. = FALSE)
  }
  init$v_sq <- as.numeric(init$v_sq)
  if (length(init$v_sq) != d) stop("init$v_sq must have length ncol(y).", call. = FALSE)
  init$v_sq <- pmax(init$v_sq, 1e-8)
  init$phi <- as.numeric(init$phi)[1]
  init$phi <- max(min(init$phi, 0.95), -0.95)
  init$p <- as.numeric(init$p)
  if (length(init$p) != J) stop("init$p must have length J.", call. = FALSE)
  init$p <- pmax(init$p, 1e-8)
  init$p <- init$p / sum(init$p)
  init
}

# ---- log-sum-exp ----
#' @keywords internal
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

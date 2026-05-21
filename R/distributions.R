# distributions.R
# Custom two-block SAD(1) distribution for paired (K = 2) data: density + RNG.
# Both blocks (the two individuals) share a single time-correlation parameter
# `phi`, while each time point keeps its own innovation variance `v_sq`
# (see paper Section 2.3). Not exported, but must exist as a package object so
# nimble can resolve it by name.

#' @keywords internal
dSADmvnorm <- nimble::nimbleFunction(
  run = function(
    x    = double(1),
    mean = double(1),
    phi  = double(0),
    v_sq = double(1),
    d1   = integer(0),
    d2   = integer(0),
    log  = integer(0, default = 0L)
  ) {
    returnType(double(0))

    x1 <- numeric(d1, init = TRUE)
    m1 <- numeric(d1, init = TRUE)
    v1 <- numeric(d1, init = TRUE)
    x2 <- numeric(d2, init = TRUE)
    m2 <- numeric(d2, init = TRUE)
    v2 <- numeric(d2, init = TRUE)

    for (jj in 1:d1) {
      x1[jj] <- x[jj]
      m1[jj] <- mean[jj]
      v1[jj] <- v_sq[jj]
    }
    for (jj in 1:d2) {
      off <- d1 + jj
      x2[jj] <- x[off]
      m2[jj] <- mean[off]
      v2[jj] <- v_sq[off]
    }

    q <- 0.0

    # Block 1 (individual 1): shared phi
    s <- 0.0
    for (jj in 1:d1) {
      s <- (x1[jj] - m1[jj]) + phi * s
      q <- q + (s * s) / v1[jj] + log(2.0 * 3.141592653589793 * v1[jj])
    }

    # Block 2 (individual 2): same shared phi
    s <- 0.0
    for (jj in 1:d2) {
      s <- (x2[jj] - m2[jj]) + phi * s
      q <- q + (s * s) / v2[jj] + log(2.0 * 3.141592653589793 * v2[jj])
    }

    logdens <- -0.5 * q
    if (log) return(logdens) else return(exp(logdens))
  }
)

#' @keywords internal
rSADmvnorm <- nimble::nimbleFunction(
  run = function(
    n    = integer(0, default = 1),
    mean = double(1),
    phi  = double(0),
    v_sq = double(1),
    d1   = integer(0),
    d2   = integer(0)
  ) {
    returnType(double(1))
    if (n != 1) print("rSADmvnorm only supports n=1")

    d <- d1 + d2
    out <- numeric(d, init = TRUE)

    s_prev <- 0.0
    for (jj in 1:d1) {
      err_j <- rnorm(1, 0, sd = sqrt(v_sq[jj]))
      x_j_minus_mean_j <- err_j - phi * s_prev
      out[jj] <- mean[jj] + x_j_minus_mean_j
      s_prev <- err_j
    }

    s_prev <- 0.0
    for (jj in 1:d2) {
      off <- d1 + jj
      err_j <- rnorm(1, 0, sd = sqrt(v_sq[off]))
      x_j_minus_mean_j <- err_j - phi * s_prev
      out[off] <- mean[off] + x_j_minus_mean_j
      s_prev <- err_j
    }

    return(out)
  }
)

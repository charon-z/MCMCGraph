# distributions.R
# 自定义 SAD(1) 两块（二元 K=2）分布：密度 + 随机数
# 不导出（internal），但必须作为包对象存在，供 nimble 通过名字找到。

#' @keywords internal
dSADmvnorm <- nimble::nimbleFunction(
  run = function(
    x    = double(1),
    mean = double(1),
    phi1 = double(0),
    phi2 = double(0),
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

    s <- 0.0
    for (jj in 1:d1) {
      s <- (x1[jj] - m1[jj]) + phi1 * s
      q <- q + (s * s) / v1[jj] + log(2.0 * 3.141592653589793 * v1[jj])
    }

    s <- 0.0
    for (jj in 1:d2) {
      s <- (x2[jj] - m2[jj]) + phi2 * s
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
    phi1 = double(0),
    phi2 = double(0),
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
      x_j_minus_mean_j <- err_j - phi1 * s_prev
      out[jj] <- mean[jj] + x_j_minus_mean_j
      s_prev <- err_j
    }

    s_prev <- 0.0
    for (jj in 1:d2) {
      off <- d1 + jj
      err_j <- rnorm(1, 0, sd = sqrt(v_sq[off]))
      x_j_minus_mean_j <- err_j - phi2 * s_prev
      out[off] <- mean[off] + x_j_minus_mean_j
      s_prev <- err_j
    }

    return(out)
  }
)

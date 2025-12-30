#' @keywords internal
.mcmcgraph_register_sad <- local({
  registered <- FALSE

  function(force = FALSE) {
    if (registered && !isTRUE(force)) return(invisible(TRUE))

    # 1) 定义 dSADmvnorm / rSADmvnorm（nimbleFunction 对象）
    #    注意：如果你已经在别的文件里定义过，就不要重复定义。
    #    最稳的做法：这里用 get0 检查，不存在才定义。
    if (is.null(get0("dSADmvnorm", envir = parent.env(environment()), inherits = FALSE))) {

      dSADmvnorm <<- nimble::nimbleFunction(
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

      rSADmvnorm <<- nimble::nimbleFunction(
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
    }

    # 2) 注册分布（nimble 需要在 build/compile 前注册）
    nimble::registerDistributions(list(
      dSADmvnorm = list(
        BUGSdist = "dSADmvnorm(mean, phi1, phi2, v_sq, d1, d2)",
        types    = c(
          "value=double(1)",
          "mean=double(1)",
          "phi1=double(0)",
          "phi2=double(0)",
          "v_sq=double(1)",
          "d1=integer(0)",
          "d2=integer(0)"
        ),
        discrete = FALSE
      )
    ))

    registered <<- TRUE
    invisible(TRUE)
  }
})

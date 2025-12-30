#' @keywords internal
.mcmcgraph_register_sad <- local({
  registered <- FALSE
  function(force = FALSE) {
    if (registered && !isTRUE(force)) return(invisible(TRUE))

    # 确保对象已经存在（否则说明 distributions.R 没被加载进来）
    if (!exists("dSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE)) {
      stop("dSADmvnorm not found in package namespace. Check R/distributions.R.")
    }
    if (!exists("rSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE)) {
      stop("rSADmvnorm not found in package namespace. Check R/distributions.R.")
    }

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

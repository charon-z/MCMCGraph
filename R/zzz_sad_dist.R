# zzz_sad_dist.R

#' @keywords internal
.mcmcgraph_expose_sad_to_global <- function() {
  # nimble 在 build/compile 时经常从 GlobalEnv 找 d/r 函数名
  if (!exists("dSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("dSADmvnorm",
           get("dSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE),
           envir = .GlobalEnv)
  }
  if (!exists("rSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("rSADmvnorm",
           get("rSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE),
           envir = .GlobalEnv)
  }
  invisible(TRUE)
}

#' Register SAD distribution in nimble (once per session)
#' @keywords internal
.mcmcgraph_register_sad <- local({
  registered <- FALSE
  function(force = FALSE) {
    if (registered && !isTRUE(force)) return(invisible(TRUE))

    if (!requireNamespace("nimble", quietly = TRUE)) {
      stop("Package 'nimble' is required. Please install it via install.packages('nimble').")
    }

    # 确保 namespace 里有 d/r
    if (!exists("dSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE))
      stop("dSADmvnorm not found in MCMCGraph namespace. Check R/distributions.R.")
    if (!exists("rSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE))
      stop("rSADmvnorm not found in MCMCGraph namespace. Check R/distributions.R.")

    # 暴露到 global 让 nimble 必然能找到
    .mcmcgraph_expose_sad_to_global()

    # ✅ nimble 1.3.0 只用它认识的字段：BUGSdist / types / discrete
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

# 兼容旧名字（如果你代码里还在调用）
#' @keywords internal
.register_sad_once <- function(...) .mcmcgraph_register_sad(...)

# zzz_sad_dist.R

#' Ensure SAD distribution functions are visible to nimble
#' @keywords internal
.mcmcgraph_expose_sad_to_global <- function() {
  # nimble 在 build/compile 时常从 GlobalEnv 找 d/r 函数名
  # 所以这里显式放一份到 global（不会锁 namespace）
  if (!exists("dSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("dSADmvnorm", get("dSADmvnorm", envir = asNamespace("MCMCGraph")), envir = .GlobalEnv)
  }
  if (!exists("rSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("rSADmvnorm", get("rSADmvnorm", envir = asNamespace("MCMCGraph")), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

#' Register SAD distribution in nimble (once per session)
#' @keywords internal
.mcmcgraph_register_sad <- local({
  registered <- FALSE
  function(force = FALSE) {
    if (registered && !isTRUE(force)) return(invisible(TRUE))

    # 1) 确认包 namespace 里确实有这俩对象
    if (!exists("dSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE)) {
      stop("dSADmvnorm not found in MCMCGraph namespace. Check R/distributions.R.")
    }
    if (!exists("rSADmvnorm", envir = asNamespace("MCMCGraph"), inherits = FALSE)) {
      stop("rSADmvnorm not found in MCMCGraph namespace. Check R/distributions.R.")
    }

    # 2) 注入到 global，避免 nimble 找不到
    .mcmcgraph_expose_sad_to_global()

    # 3) 注册分布（明确告诉 nimble：dfunc/rfunc 分别是谁）
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
        discrete = FALSE,
        pqAvail  = FALSE,
        # 关键：把函数名显式写进去
        dfunc = "dSADmvnorm",
        rfunc = "rSADmvnorm"
      )
    ))

    registered <<- TRUE
    invisible(TRUE)
  }
})

# 兼容你之前写过的名字（如果 run_mcmc_graph.R 里还在调用旧名）
#' @keywords internal
.register_sad_once <- function(...) .mcmcgraph_register_sad(...)

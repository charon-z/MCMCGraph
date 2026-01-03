# zzz_sad_dist.R
# 状态环境（避免 namespace 锁定绑定问题）
.mcmcgraph_state <- new.env(parent = emptyenv())
.mcmcgraph_state$registered <- FALSE

# 把分布函数暴露到 GlobalEnv（nimble build/compile 更稳）
.mcmcgraph_expose_sad_to_global <- function(force = FALSE) {
  if (force || !exists("dSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("dSADmvnorm", get("dSADmvnorm", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  if (force || !exists("rSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("rSADmvnorm", get("rSADmvnorm", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

# 注册（只需一次）
.mcmcgraph_register_sad <- function(force = FALSE) {
  # 1) 保证 d/r 在 global（nimble 查找用）
  .mcmcgraph_expose_sad_to_global(force = force)

  # 2) 只注册一次（除非 force）
  if (!force && isTRUE(.mcmcgraph_state$registered)) return(invisible(TRUE))

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

  .mcmcgraph_state$registered <- TRUE
  invisible(TRUE)
}

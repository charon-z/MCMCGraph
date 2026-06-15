# zzz_sad_dist.R
# State environment (avoids namespace locked-binding issues at registration time).
.mcmcgraph_state <- new.env(parent = emptyenv())
.mcmcgraph_state$registered <- FALSE

# Expose the distribution functions into GlobalEnv (more robust for nimble
# build/compile, which resolves the d/r functions by name).
.mcmcgraph_expose_sad_to_global <- function(force = FALSE) {
  if (force || !exists("dSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("dSADmvnorm", get("dSADmvnorm", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  if (force || !exists("rSADmvnorm", envir = .GlobalEnv, inherits = FALSE)) {
    assign("rSADmvnorm", get("rSADmvnorm", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  if (force || !exists("dSADmvnorm2", envir = .GlobalEnv, inherits = FALSE)) {
    assign("dSADmvnorm2", get("dSADmvnorm2", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  if (force || !exists("rSADmvnorm2", envir = .GlobalEnv, inherits = FALSE)) {
    assign("rSADmvnorm2", get("rSADmvnorm2", envir = parent.env(environment())), envir = .GlobalEnv)
  }
  invisible(TRUE)
}

# Register the SAD distribution with nimble (only needs to happen once).
.mcmcgraph_register_sad <- function(force = FALSE) {
  # 1) make sure d/r live in GlobalEnv (nimble looks them up by name)
  .mcmcgraph_expose_sad_to_global(force = force)

  # 2) register only once (unless force = TRUE)
  if (!force && isTRUE(.mcmcgraph_state$registered)) return(invisible(TRUE))

  nimble::registerDistributions(list(
    dSADmvnorm = list(
      BUGSdist = "dSADmvnorm(mean, phi, v_sq, d1, d2)",
      types    = c(
        "value=double(1)",
        "mean=double(1)",
        "phi=double(0)",
        "v_sq=double(1)",
        "d1=integer(0)",
        "d2=integer(0)"
      ),
      discrete = FALSE
    ),
    dSADmvnorm2 = list(
      BUGSdist = "dSADmvnorm2(mean, phi, psi, v_sq, d1, d2)",
      types    = c(
        "value=double(1)",
        "mean=double(1)",
        "phi=double(0)",
        "psi=double(0)",
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

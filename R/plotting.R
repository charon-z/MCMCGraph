#' Plot trace + density for selected parameters
#'
#' @param result mcmcgraph_result from run_mcmc_binary()
#' @param params character vector of parameter names (regex allowed)
#' @param max_params max number of params to plot
#' @export
plot_trace_density <- function(result, params = c("^phi", "^p\\["),
                               max_params = 6) {
  if (!inherits(result, "mcmcgraph_result")) stop("result must be from run_mcmc_binary().")
  samp <- result$posterior_samples$all_params
  if (is.null(samp) || nrow(samp) < 5) stop("No samples found.")

  cn <- colnames(samp)
  keep <- character(0)
  for (p in params) keep <- c(keep, grep(p, cn, value = TRUE))
  keep <- unique(keep)
  if (length(keep) == 0) stop("No parameters matched.")
  keep <- keep[1:min(length(keep), max_params)]

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(length(keep), 2), mar = c(3, 3, 2, 1))

  for (nm in keep) {
    x <- samp[, nm]
    plot(x, type = "l", xlab = "iter", ylab = nm, main = paste0("Trace: ", nm))
    plot(density(x), xlab = nm, main = paste0("Density: ", nm))
  }
  invisible(keep)
}

#' Collect evaluation files produced by user scripts (optional helper)
#'
#' @param eval_list list of evaluation objects OR data.frame
#' @return data.frame with J and BIC
#' @export
as_eval_df <- function(eval_list) {
  if (is.data.frame(eval_list)) return(eval_list)
  if (!is.list(eval_list)) stop("eval_list must be a list or data.frame.")
  rows <- lapply(eval_list, function(e) {
    if (is.null(e$J) || is.null(e$BIC)) stop("Each eval object must contain $J and $BIC.")
    data.frame(J = e$J, BIC = e$BIC, loglik = if (!is.null(e$loglik)) e$loglik else NA_real_)
  })
  do.call(rbind, rows)
}

#' Plot BIC vs J
#'
#' @param eval_df data.frame with columns J and BIC
#' @export
plot_bic <- function(eval_df) {
  df <- as_eval_df(eval_df)
  df <- df[order(df$J), , drop = FALSE]
  plot(df$J, df$BIC, type = "b", xlab = "J", ylab = "BIC", main = "BIC vs J")
  invisible(df)
}

#' Plot mixing-proportion MCMC traces
#'
#' Draw one trace panel for each selected mixture component. The ordinate is
#' labelled explicitly as the mixing proportion, while individual panels are
#' identified by the corresponding mathematical symbol, \ifelse{html}{\out{
#' &pi;<sub>j</sub>}}{\eqn{\pi_j}}. By default, panels use separate y-axis
#' ranges so that convergence remains visible for both large and small
#' components.
#'
#' @param result A `mcmcgraph_result` returned by [run_mcmc_binary()].
#' @param components Integer vector selecting mixture components. `NULL` plots
#'   all components.
#' @param ncol Number of panel columns.
#' @param free_y If `TRUE`, use a separate y-axis range for each component. If
#'   `FALSE`, use one common range across all selected components.
#' @param show_running_mean Draw the cumulative running mean.
#' @param show_posterior_mean Draw a dashed horizontal posterior-mean line.
#' @param show_legend Add a compact legend to the first panel.
#' @param trace_col,runmean_col,mean_col Colours for posterior draws, running
#'   means and posterior means.
#' @return Invisibly, a data frame containing the selected component numbers
#'   and posterior means. Called primarily for its plotting side effect.
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' fit <- run_mcmc_binary(example_binary$y, J = 3,
#'                        times = example_binary$times,
#'                        niter = 300, seed = 1)
#' plot_mixing_trace(fit)
#' }
plot_mixing_trace <- function(
    result,
    components = NULL,
    ncol = 3L,
    free_y = TRUE,
    show_running_mean = TRUE,
    show_posterior_mean = TRUE,
    show_legend = TRUE,
    trace_col = "grey70",
    runmean_col = "#D55E00",
    mean_col = "#0072B2"
) {
  if (!inherits(result, "mcmcgraph_result")) {
    stop("result must be from run_mcmc_binary().")
  }

  p <- result$posterior_samples$p
  if (is.null(p)) stop("No mixing-proportion samples found in result.")
  p <- as.matrix(p)
  if (nrow(p) < 2L || ncol(p) < 1L) {
    stop("Mixing-proportion samples must contain at least two draws.")
  }

  J <- ncol(p)
  if (is.null(components)) components <- seq_len(J)
  if (!is.numeric(components) || anyNA(components) ||
      any(components != as.integer(components)) ||
      any(components < 1L | components > J)) {
    stop("components must contain valid integer component numbers.")
  }
  components <- unique(as.integer(components))
  if (!length(components)) stop("At least one component must be selected.")

  ncol <- as.integer(ncol)
  if (length(ncol) != 1L || is.na(ncol) || ncol < 1L) {
    stop("ncol must be a positive integer.")
  }

  selected <- p[, components, drop = FALSE]
  if (any(!is.finite(selected))) {
    stop("Mixing-proportion samples must be finite.")
  }

  model_info <- result$model_info
  thin <- if (!is.null(model_info$thin)) as.integer(model_info$thin) else 1L
  burnin_fraction <- if (!is.null(model_info$burnin_frac)) {
    as.numeric(model_info$burnin_frac)
  } else {
    0
  }
  total_iterations <- if (!is.null(model_info$niter)) {
    as.integer(model_info$niter)
  } else {
    NA_integer_
  }
  if (is.finite(total_iterations) && is.finite(burnin_fraction)) {
    first_iteration <- floor(total_iterations * burnin_fraction) + thin
    iteration <- first_iteration + (seq_len(nrow(selected)) - 1L) * thin
  } else {
    iteration <- seq_len(nrow(selected))
  }

  common_ylim <- range(selected)
  panel_rows <- ceiling(length(components) / ncol)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(panel_rows, ncol),
    mar = c(2.8, 3.8, 2.2, 1.0),
    oma = c(3.2, 4.8, 0.6, 0.4)
  )

  for (panel_index in seq_along(components)) {
    component <- components[panel_index]
    values <- selected[, panel_index]
    panel_ylim <- if (isTRUE(free_y)) range(values) else common_ylim
    if (diff(panel_ylim) == 0) {
      padding <- max(abs(panel_ylim[1L]) * 0.04, 1e-6)
      panel_ylim <- panel_ylim + c(-padding, padding)
    }

    graphics::plot(
      iteration, values,
      type = "l", col = trace_col, lwd = 0.7,
      axes = FALSE, xlab = "", ylab = "", ylim = panel_ylim
    )
    graphics::axis(1)
    graphics::axis(2, las = 1)
    graphics::box(col = "grey35")
    graphics::title(
      main = bquote(pi[.(component)]),
      adj = 0, line = 0.35, font.main = 1
    )

    if (isTRUE(show_running_mean)) {
      graphics::lines(
        iteration, cumsum(values) / seq_along(values),
        col = runmean_col, lwd = 1.8
      )
    }
    if (isTRUE(show_posterior_mean)) {
      graphics::abline(h = mean(values), col = mean_col, lwd = 1.3, lty = 2)
    }

    if (panel_index == 1L && isTRUE(show_legend)) {
      legend_labels <- "Posterior draw"
      legend_colours <- trace_col
      legend_types <- 1
      legend_widths <- 0.7
      if (isTRUE(show_running_mean)) {
        legend_labels <- c(legend_labels, "Running mean")
        legend_colours <- c(legend_colours, runmean_col)
        legend_types <- c(legend_types, 1)
        legend_widths <- c(legend_widths, 1.8)
      }
      if (isTRUE(show_posterior_mean)) {
        legend_labels <- c(legend_labels, "Posterior mean")
        legend_colours <- c(legend_colours, mean_col)
        legend_types <- c(legend_types, 2)
        legend_widths <- c(legend_widths, 1.3)
      }
      graphics::legend(
        "topright", legend = legend_labels, col = legend_colours,
        lty = legend_types, lwd = legend_widths, bty = "n", cex = 0.72
      )
    }
  }

  empty_panels <- panel_rows * ncol - length(components)
  if (empty_panels > 0L) {
    for (unused in seq_len(empty_panels)) graphics::plot.new()
  }

  graphics::mtext(
    "Post-burn-in MCMC iteration", side = 1, outer = TRUE, line = 1.5
  )
  graphics::mtext(
    expression("Mixing proportion (" * pi[j] * ")"),
    side = 2, outer = TRUE, line = 2.6
  )

  invisible(data.frame(
    component = components,
    posterior_mean = colMeans(selected),
    row.names = NULL
  ))
}

#' Plot trace + density for selected parameters
#'
#' @param result mcmcgraph_result from run_mcmc_binary()
#' @param params character vector of parameter names (regex allowed)
#' @param max_params max number of params to plot
#' @return (Invisibly) the character vector of parameter names that were
#'   plotted. Called for its side effect of drawing trace and density panels.
#' @export
#' @examples
#' \donttest{
#' data(example_binary)
#' fit <- run_mcmc_binary(example_binary$y, J = 3, times = example_binary$times,
#'                        niter = 300, seed = 1)
#' plot_trace_density(fit, params = c("^phi"))
#' }
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
#' @examples
#' as_eval_df(data.frame(J = 2:4, BIC = c(120, 100, 110), loglik = c(-50, -40, -42)))
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
#' @return (Invisibly) the ordered data.frame that was plotted. Called for its
#'   side effect of drawing the BIC-versus-J curve.
#' @export
#' @examples
#' plot_bic(data.frame(J = 2:5, BIC = c(140, 110, 115, 130)))
plot_bic <- function(eval_df) {
  df <- as_eval_df(eval_df)
  df <- df[order(df$J), , drop = FALSE]
  plot(df$J, df$BIC, type = "b", xlab = "J", ylab = "BIC", main = "BIC vs J")
  invisible(df)
}

make_mock_mcmc_result <- function() {
  set.seed(42)
  p <- matrix(stats::runif(120), nrow = 40, ncol = 3)
  p <- p / rowSums(p)
  colnames(p) <- paste0("p[", seq_len(ncol(p)), "]")
  structure(
    list(
      posterior_samples = list(p = p),
      model_info = list(niter = 80L, thin = 1L, burnin_frac = 0.5)
    ),
    class = "mcmcgraph_result"
  )
}

test_that("plot_mixing_trace returns component summaries", {
  result <- make_mock_mcmc_result()
  output <- tempfile(fileext = ".png")
  grDevices::png(output, width = 900, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)

  summary <- plot_mixing_trace(
    result, components = c(1, 3), ncol = 2, show_legend = FALSE
  )

  expect_s3_class(summary, "data.frame")
  expect_equal(summary$component, c(1L, 3L))
  expect_equal(
    summary$posterior_mean,
    unname(colMeans(result$posterior_samples$p[, c(1, 3)]))
  )
})

test_that("plot_mixing_trace validates inputs", {
  result <- make_mock_mcmc_result()
  expect_error(plot_mixing_trace(list()), "run_mcmc_binary")
  expect_error(plot_mixing_trace(result, components = 0), "valid integer")
  expect_error(plot_mixing_trace(result, components = 4), "valid integer")
  expect_error(plot_mixing_trace(result, ncol = 0), "positive integer")

  result$posterior_samples$p <- NULL
  expect_error(plot_mixing_trace(result), "No mixing-proportion samples")
})

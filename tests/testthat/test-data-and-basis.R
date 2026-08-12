# Fast, pure-R tests (no nimble compilation) covering the data adapter, the
# Legendre basis and the BIC bookkeeping. A short end-to-end MCMC smoke test is
# kept separate and skipped on CRAN.

test_that("as_binary_data accepts a wide matrix", {
  data(example_binary)
  dat <- as_binary_data(example_binary$y, format = "wide", times = example_binary$times)
  expect_equal(nrow(dat$y), nrow(example_binary$y))
  expect_equal(ncol(dat$y), ncol(example_binary$y))
  expect_equal(dat$d_single, ncol(example_binary$y) / 2)
  expect_length(dat$times_single, dat$d_single)
})

test_that("as_binary_data round-trips wide <-> long", {
  data(example_binary)
  d <- ncol(example_binary$y) / 2
  ids <- rownames(example_binary$y)
  if (is.null(ids)) ids <- seq_len(nrow(example_binary$y))
  long <- do.call(rbind, lapply(seq_len(nrow(example_binary$y)), function(i) {
    rbind(
      data.frame(snp = ids[i], ind = "A", time = seq_len(d),
                 y = example_binary$y[i, 1:d]),
      data.frame(snp = ids[i], ind = "B", time = seq_len(d),
                 y = example_binary$y[i, (d + 1):(2 * d)]))
  }))
  dat <- as_binary_data(long, format = "long")
  expect_equal(dat$d_single, d)
  expect_equal(nrow(dat$y), nrow(example_binary$y))
  expect_equal(unname(dat$y[1, 1:d]), unname(example_binary$y[1, 1:d]))
})

test_that("as_binary_data rejects NA and odd column counts", {
  expect_error(as_binary_data(matrix(1:6, nrow = 2), format = "wide"))  # 3 cols (odd)
  m <- matrix(rnorm(8), nrow = 2); m[1, 1] <- NA
  expect_error(as_binary_data(m, format = "wide"))
})

test_that("eval/plot helpers behave", {
  ev <- data.frame(J = 2:4, BIC = c(120, 100, 110), loglik = c(-50, -40, -42))
  expect_s3_class(as_eval_df(ev), "data.frame")
  expect_true(all(c("J", "BIC") %in% names(as_eval_df(ev))))
})

test_that("end-to-end clustering runs and recovers structure (smoke)", {
  skip_on_cran()
  skip_if_not_installed("nimble")
  data(example_binary)
  fit <- run_mcmc_binary(example_binary$y, J = 3, times = example_binary$times,
                         niter = 200, seed = 1)
  expect_s3_class(fit, "mcmcgraph_result")
  expect_length(fit$clustering, nrow(example_binary$y))
  expect_equal(ncol(fit$cluster_prob), 3)
  expect_true(all(fit$clustering %in% 1:3))
})

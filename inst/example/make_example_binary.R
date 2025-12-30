# inst/example/make_example_binary.R

set.seed(20251230)

# ---- knobs ----
n <- 100          # 小一点，跑示例很快
d <- 16          # 示例固定16（你的方法支持任意d，示例选最常见）
times <- 1:d

# LOP order fixed = 4 -> q = 5 -> P = 10
order <- 4
q <- order + 1

# ---- build Z0 (same as package logic, copy here to generate data) ----
calc_Z_legendre <- function(times, order = 4) {
  d <- length(times)
  q <- order + 1L
  Z <- matrix(0, d, q)
  tmin <- min(times); tmax <- max(times)
  tx <- if (tmax == tmin) rep(0, d) else (-1 + 2 * (times - tmin) / (tmax - tmin))
  if (order > 4) stop("order > 4 not supported in this template.")
  for (i in seq_len(d)) {
    x <- tx[i]
    vals <- c(1,
              x,
              0.5 * (3*x^2 - 1),
              0.5 * (5*x^3 - 3*x),
              0.125 * (35*x^4 - 30*x^2 + 3))
    Z[i, ] <- vals[1:q]
  }
  Z
}

make_Z0_binary <- function(times_single, order = 4) {
  Z1 <- calc_Z_legendre(times_single, order = order)  # d x q
  d <- nrow(Z1); q <- ncol(Z1)
  Zb <- matrix(0, nrow = 2*d, ncol = 2*q)
  Zb[1:d, 1:q] <- Z1
  Zb[(d+1):(2*d), (q+1):(2*q)] <- Z1
  Zb
}

Z0_binary <- make_Z0_binary(times, order = order) # (2d) x 10

# ---- simulate "clustered" mean curves (optional, just for realism) ----
J_true <- 3
z_true <- sample(1:J_true, n, replace = TRUE, prob = c(0.4, 0.35, 0.25))

# beta_true: J x 10
beta_true <- matrix(0, nrow = J_true, ncol = 2*q)
for (j in 1:J_true) {
  # individual 1 coefficients
  beta_true[j, 1:q] <- rnorm(q, mean = (j-2)*0.4, sd = 0.4)
  # individual 2 coefficients (slightly shifted)
  beta_true[j, (q+1):(2*q)] <- rnorm(q, mean = (2-j)*0.3, sd = 0.4)
}

mu_mat <- beta_true %*% t(Z0_binary)  # J x (2d)

# ---- simulate SAD-like noise by two AR(1) recursions (phi1/phi2) ----
phi1 <- 0.35
phi2 <- 0.20
v_sq <- rep(0.6, 2*d)  # constant variance for simplicity

sim_two_block <- function(mu, phi1, phi2, v_sq, d) {
  stopifnot(length(mu) == 2*d, length(v_sq) == 2*d)
  y <- numeric(2*d)

  # block1
  s <- 0
  for (t in 1:d) {
    e <- rnorm(1, 0, sd = sqrt(v_sq[t]))
    # reverse of your likelihood recursion; just a reasonable correlated noise
    s <- e - phi1 * s
    y[t] <- mu[t] + s
  }
  # block2
  s <- 0
  for (t in 1:d) {
    idx <- d + t
    e <- rnorm(1, 0, sd = sqrt(v_sq[idx]))
    s <- e - phi2 * s
    y[idx] <- mu[idx] + s
  }
  y
}

Y <- matrix(0, nrow = n, ncol = 2*d)
for (i in 1:n) {
  mu_i <- mu_mat[z_true[i], ]
  Y[i, ] <- sim_two_block(mu_i, phi1, phi2, v_sq, d)
}

# ---- add informative column names ----
cn1 <- sprintf("i1_t%02d", 1:d)
cn2 <- sprintf("i2_t%02d", 1:d)
colnames(Y) <- c(cn1, cn2)
rownames(Y) <- sprintf("snp_%03d", 1:n)

example_binary <- list(
  y = Y,
  times = times,
  d_single = d,
  order = order,
  truth = list(J_true = J_true, z_true = z_true)
)

# ---- write a human-readable template/example file (optional) ----
dir.create("inst/example", recursive = TRUE, showWarnings = FALSE)
write.table(
  example_binary$y,
  file = "inst/example/example_binary_small.tsv",
  sep = "\t", quote = FALSE, col.names = NA
)

# ---- save into R package data/ as .rda ----
# Run this script from package root so usethis can find package.
if (!requireNamespace("usethis", quietly = TRUE)) install.packages("usethis")
usethis::use_data(example_binary, overwrite = TRUE)

message("Done. Saved:\n- inst/example/example_binary_small.tsv\n- data/example_binary.rda (load via data('example_binary'))")

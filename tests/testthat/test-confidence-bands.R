# Tests for confidence_bands() (bd-3c0.13, bd-2qe.13)

make_fitted_mean <- function() {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  estimate(m)
}

make_fitted_mean_variance <- function() {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y
  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  estimate(m)
}

# confidence_bands() via MEstimator -------------------------------------------

test_that("confidence_bands() returns p-by-2 matrix", {
  m <- make_fitted_mean_variance()
  cb <- confidence_bands(m, seed = 42)
  expect_true(is.matrix(cb))
  expect_equal(nrow(cb), length(m@theta))
  expect_equal(ncol(cb), 2)
  expect_true(all(c("lower", "upper") %in% colnames(cb)))
})

test_that("confidence_bands() for single parameter", {
  m <- make_fitted_mean()
  cb <- confidence_bands(m, seed = 42)
  expect_equal(nrow(cb), 1)
  expect_equal(ncol(cb), 2)
})

test_that("confidence_bands() are wider than confidence intervals", {
  m <- make_fitted_mean_variance()
  ci <- confidence_intervals(m, alpha = 0.05)
  cb <- confidence_bands(m, alpha = 0.05, seed = 42)
  # Bands should be at least as wide as intervals for each parameter
  ci_widths <- ci[, "upper"] - ci[, "lower"]
  cb_widths <- cb[, "upper"] - cb[, "lower"]
  expect_true(all(cb_widths >= ci_widths - 1e-10))
})

test_that("confidence_bands() with bonferroni method works", {
  m <- make_fitted_mean_variance()
  cb <- confidence_bands(m, method = "bonferroni")
  expect_true(is.matrix(cb))
  expect_equal(nrow(cb), 2)
  expect_equal(ncol(cb), 2)
})

test_that("confidence_bands() bonferroni wider than pointwise CI", {
  m <- make_fitted_mean_variance()
  ci <- confidence_intervals(m, alpha = 0.05)
  cb <- confidence_bands(m, method = "bonferroni", alpha = 0.05)
  ci_widths <- ci[, "upper"] - ci[, "lower"]
  cb_widths <- cb[, "upper"] - cb[, "lower"]
  expect_true(all(cb_widths >= ci_widths - 1e-10))
})

test_that("confidence_bands() with subset parameter", {
  m <- make_fitted_mean_variance()
  cb <- confidence_bands(m, subset = 1L, seed = 42)
  expect_equal(nrow(cb), 1)
  expect_equal(ncol(cb), 2)
})

test_that("confidence_bands() errors before estimation", {
  y <- c(1, 2, 3)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(confidence_bands(m))
})

test_that("confidence_bands() with seed is reproducible", {
  m <- make_fitted_mean_variance()
  cb1 <- confidence_bands(m, seed = 123)
  cb2 <- confidence_bands(m, seed = 123)
  expect_equal(cb1, cb2)
})

test_that("confidence_bands() with different seeds gives different results", {
  m <- make_fitted_mean_variance()
  cb1 <- confidence_bands(m, seed = 1)
  cb2 <- confidence_bands(m, seed = 2)
  # Very unlikely to be exactly equal with different seeds
  expect_false(identical(cb1, cb2))
})

# n_draws default and forwarding ---------------------------------------------

test_that("confidence_bands() default n_draws is 1e5", {
  # deli deliberately keeps 1e5 rather than Python's 1e6 estimator default.
  expect_identical(formals(args(confidence_bands))$n_draws, 100000L)
  expect_identical(formals(compute_confidence_bands)$n_draws, 100000L)
})

test_that("confidence_bands() forwards the default n_draws to the sup-t draw", {
  m <- make_fitted_mean_variance()
  # With a shared seed, the default path must reproduce an explicit 1e5 call
  # exactly (proving 1e5 is the default and that it is forwarded)...
  default_cb <- confidence_bands(m, seed = 7)
  explicit_1e5 <- compute_confidence_bands(
    m@theta,
    m@variance,
    seed = 7,
    n_draws = 100000L
  )
  expect_equal(default_cb, explicit_1e5)
  # ...and must differ from a 1e6 call, confirming the default is not 1e6.
  explicit_1e6 <- compute_confidence_bands(
    m@theta,
    m@variance,
    seed = 7,
    n_draws = 1000000L
  )
  expect_false(isTRUE(all.equal(default_cb, explicit_1e6)))
})

test_that("confidence_bands() forwards an explicit n_draws", {
  m <- make_fitted_mean_variance()
  forwarded <- confidence_bands(m, seed = 11, n_draws = 1000000L)
  direct <- compute_confidence_bands(
    m@theta,
    m@variance,
    seed = 11,
    n_draws = 1000000L
  )
  expect_equal(forwarded, direct)
})

# compute_confidence_bands() standalone function ------------------------------

test_that("compute_confidence_bands() works with raw theta and covariance", {
  theta <- c(1, 2)
  covariance <- matrix(c(0.1, 0.02, 0.02, 0.2), nrow = 2)
  cb <- compute_confidence_bands(theta, covariance, seed = 42)
  expect_true(is.matrix(cb))
  expect_equal(nrow(cb), 2)
  expect_equal(ncol(cb), 2)
})

test_that("bonferroni critical value matches qnorm(1 - alpha/(2k))", {
  theta <- c(1, 2, 3)
  covariance <- diag(c(0.1, 0.2, 0.3))
  cb <- compute_confidence_bands(theta, covariance, method = "bonferroni")
  se <- sqrt(diag(covariance))
  # Bonferroni CV = qnorm(1 - alpha / (2*k))
  k <- length(theta)
  cv_expected <- qnorm(1 - 0.05 / (2 * k))
  expected_lower <- theta - cv_expected * se
  expected_upper <- theta + cv_expected * se
  expect_equal(cb[, "lower"], expected_lower, tolerance = 1e-10)
  expect_equal(cb[, "upper"], expected_upper, tolerance = 1e-10)
})

test_that("confidence_bands() narrower alpha means wider bands", {
  m <- make_fitted_mean_variance()
  cb_95 <- confidence_bands(m, alpha = 0.05, method = "bonferroni")
  cb_99 <- confidence_bands(m, alpha = 0.01, method = "bonferroni")
  width_95 <- cb_95[, "upper"] - cb_95[, "lower"]
  width_99 <- cb_99[, "upper"] - cb_99[, "lower"]
  expect_true(all(width_99 > width_95))
})

# compute_confidence_bands() input validation ---------------------------------
#
# Python's confidence_bands raises ValueError when the parameter vector and the
# covariance matrix disagree on dimension, and when alpha is not strictly between
# 0 and 1. Without these checks, R silently recycles a mismatched theta over the
# covariance (a length-4 theta on a 2x2 covariance returns four rows of garbage),
# warns-but-returns recycled values for a length-3 theta, returns (-Inf, Inf)
# bands at alpha = 0, and returns NaN bands at alpha = -0.5. These tests pin
# informative aborts so a bad shape or alpha fails loudly.

test_that("compute_confidence_bands() aborts when theta is longer than covariance", {
  # A length-4 theta over a 2x2 covariance silently recycles to four rows.
  expect_error(
    compute_confidence_bands(
      c(1, 2, 3, 4),
      diag(2),
      method = "bonferroni"
    ),
    regexp = "dimension"
  )
})

test_that("compute_confidence_bands() aborts when theta is shorter than covariance", {
  expect_error(
    compute_confidence_bands(
      c(1, 2, 3),
      diag(2),
      method = "bonferroni"
    ),
    regexp = "dimension"
  )
})

test_that("compute_confidence_bands() aborts on a non-square covariance matrix", {
  expect_error(
    compute_confidence_bands(
      c(1, 2),
      matrix(seq_len(6), nrow = 2, ncol = 3),
      method = "bonferroni"
    ),
    regexp = "square"
  )
})

test_that("compute_confidence_bands() aborts on alpha = 0", {
  expect_error(
    compute_confidence_bands(
      c(1, 2),
      diag(2),
      alpha = 0,
      method = "bonferroni"
    ),
    regexp = "between 0 and 1"
  )
})

test_that("compute_confidence_bands() aborts on a negative alpha", {
  expect_error(
    compute_confidence_bands(
      c(1, 2),
      diag(2),
      alpha = -0.5,
      method = "bonferroni"
    ),
    regexp = "between 0 and 1"
  )
})

test_that("compute_confidence_bands() aborts on alpha = 1", {
  expect_error(
    compute_confidence_bands(
      c(1, 2),
      diag(2),
      alpha = 1,
      method = "bonferroni"
    ),
    regexp = "between 0 and 1"
  )
})

# numeric confidence_bands() subset -------------------------------------------
#
# The numeric method accepts subset in its signature but never applies it, so
# subset = c(1, 2) on a 3-parameter input returns all three rows. The estimator
# method applies subset correctly. These tests pin the numeric method to the
# same behavior: subset the parameter vector and covariance before computing.

test_that("numeric confidence_bands() honors the subset argument", {
  theta <- c(1, 2, 3)
  covariance <- diag(c(0.1, 0.2, 0.3))
  cb <- confidence_bands(
    theta,
    covariance = covariance,
    subset = c(1, 2),
    method = "bonferroni"
  )
  expect_equal(nrow(cb), 2)
  expect_equal(ncol(cb), 2)
})

test_that("numeric confidence_bands() subset matches bands on the reduced inputs", {
  theta <- c(1, 2, 3)
  covariance <- diag(c(0.1, 0.2, 0.3))
  subset <- c(1, 3)
  cb <- confidence_bands(
    theta,
    covariance = covariance,
    subset = subset,
    method = "bonferroni"
  )
  reduced <- compute_confidence_bands(
    theta[subset],
    covariance[subset, subset, drop = FALSE],
    method = "bonferroni"
  )
  expect_equal(cb, reduced)
})

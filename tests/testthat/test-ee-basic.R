# Tests for basic estimating equations using the MEstimator pipeline
# Each test defines a psi function inline, creates an MEstimator, calls
# estimate(), and validates theta and variance against Python fixtures.

test_that("ee_mean: mean of y matches Python", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  # Estimating equation for the mean: psi(theta) = y_i - theta

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  expect_equal(unname(m@variance), ref$variance, tolerance = 1e-4)
})

test_that("ee_mean_variance: mean and variance match Python", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  # Stacked estimating equations: mean and variance
  psi <- function(theta) {
    rbind(
      y - theta[1], # EE for the mean
      (y - theta[1])^2 - theta[2] # EE for the variance
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  # Compare diagonal of variance matrix
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-4)
})

test_that("ee_mean_robust with Huber loss matches Python", {
  ref <- load_fixture("ee_mean_robust_huber")
  y <- ref$y
  k <- ref$k

  # Robust mean using Huber loss function
  psi <- function(theta) {
    matrix(robust_loss_functions(y - theta[1], "huber", k = k), nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  expect_equal(unname(m@variance), ref$variance, tolerance = 1e-4)
})

test_that("ee_mean_geometric with log_theta = FALSE matches the analytic log geometric mean", {
  # Python cross-validation of ee_mean_geometric lives in test-python-fixtures.R
  # (the default log_theta = TRUE parameterization). Here we check the
  # log_theta = FALSE parameterization, where theta is the log of the geometric
  # mean, against its analytic value mean(log(y)).
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  psi <- function(theta) {
    ee_mean_geometric(theta, y = y, log_theta = FALSE)
  }

  m <- MEstimator(stacked_equations = psi, init = c(1))
  m <- estimate(m)

  # The log geometric mean of 1:9 is mean(log(1:9))
  expected_log_gm <- mean(log(y))
  expect_equal(unname(m@theta[1]), expected_log_gm, tolerance = 1e-4)
})

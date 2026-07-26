# Tests for ee_ridge_regression() estimating equations using Python fixtures
# ee_ridge_regression() implements L2-penalized regression. It takes theta,
# X, y, model, penalty, and optional weights/center/offset arguments.

test_that("ee_ridge_regression matches Python Delicatessen theta and variance", {
  ref <- load_fixture("ee_ridge_regression")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty

  # Ridge regression via ee_ridge_regression wrapper
  psi <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty
    )
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_ridge_regression pipeline results match full fixture", {
  ref <- load_fixture("ee_ridge_regression")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty

  psi <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty
    )
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Check asymptotic variance
  expect_equal(
    unname(diag(m@asymptotic_variance)),
    diag(ref$asymptotic_variance),
    tolerance = 1e-3
  )

  # Check bread matrix
  expect_equal(unname(m@bread), ref$bread, tolerance = 1e-3)

  # Check confidence intervals
  ci <- confidence_intervals(m)
  expect_equal(unname(ci[, "lower"]), ref$ci_lower, tolerance = 1e-3)
  expect_equal(unname(ci[, "upper"]), ref$ci_upper, tolerance = 1e-3)
})

test_that("ee_ridge_regression with zero penalty matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Ridge regression with penalty = 0 should be identical to unpenalized
  # linear regression
  psi_ridge <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = 0
    )
  }

  psi_linear <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }

  m_ridge <- MEstimator(stacked_equations = psi_ridge, init = init)
  m_ridge <- estimate(m_ridge)

  m_linear <- MEstimator(stacked_equations = psi_linear, init = init)
  m_linear <- estimate(m_linear)

  expect_equal(m_ridge@theta, m_linear@theta, tolerance = 1e-6)
  expect_equal(
    diag(m_ridge@variance),
    diag(m_linear@variance),
    tolerance = 1e-6
  )
})

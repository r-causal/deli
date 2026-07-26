# Tests for regression estimating equations pipeline using Python fixtures

test_that("ee_regression linear matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_linear")

  # fromJSON with simplifyVector produces a matrix for X directly
  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Linear regression estimating equation: psi(theta) = X_i * (y_i - X_i^T theta)
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - Xb
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_regression logistic matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Inverse logit (expit) link function
  inverse_logit <- function(x) {
    1 / (1 + exp(-x))
  }

  # Logistic regression estimating equation
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - inverse_logit(Xb)
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_regression poisson matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  n <- length(y)

  # Poisson regression estimating equation (exp link)
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - exp(Xb)
    t(X * as.numeric(residuals)) # p-by-n
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

test_that("ee_ridge_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_ridge_regression")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty
  n <- length(y)

  # Ridge regression: linear regression + L2 penalty
  # Python: penalty_terms = (penalty / n) * theta (for gamma=2, center=0)
  # Return: ((y - pred_y) * X).T - penalty_terms[:, None]
  psi <- function(theta) {
    Xb <- X %*% theta
    residuals <- y - Xb

    # Regression score: p-by-n matrix
    score <- t(X * as.numeric(residuals))

    # L2 penalty term: (penalty / n) * theta, broadcast across all n columns
    penalty_terms <- (penalty / n) * theta
    score - penalty_terms # recycled across columns
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  ref_var <- ref$variance
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref_var), tolerance = 1e-3)
})

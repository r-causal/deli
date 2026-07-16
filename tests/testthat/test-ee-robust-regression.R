# Tests for ee_robust_regression() estimating equations using Python fixtures
# ee_robust_regression() implements robust regression with various loss
# functions (Huber, Tukey's biweight, Andrew's Sine, etc.). It takes theta,
# X, y, model, k (tuning parameter), and loss arguments.

test_that("ee_robust_regression Huber loss matches Python Delicatessen", {
  ref <- load_fixture("ee_robust_regression_huber")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  k <- ref$k

  # Huber robust regression: applies Huber loss function to residuals
  # f_k(r) = r if |r| < k, sign(r) * k otherwise
  # Estimating equation: f_k(y - X %*% theta) * X
  psi <- function(theta) {
    ee_robust_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      k = k,
      loss = "huber"
    )
  }

  # Huber clips residuals beyond k, so from a starting point where every
  # residual is clipped the Jacobian is singular and Newton solvers stall.
  # Python delicatessen solves with scipy's Levenberg-Marquardt; the "lm"
  # solver wraps the same MINPACK routine and converges from the origin.
  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m, solver = "lm")

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-5)
})

test_that("ee_robust_regression Huber loss pipeline results match full fixture", {
  ref <- load_fixture("ee_robust_regression_huber")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  k <- ref$k

  psi <- function(theta) {
    ee_robust_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      k = k,
      loss = "huber"
    )
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m, solver = "lm")

  # Check asymptotic variance
  ref_avar <- ref$asymptotic_variance
  expect_equal(unname(diag(m@asymptotic_variance)), diag(ref_avar), tolerance = 1e-5)

  # Check bread matrix
  ref_bread <- ref$bread
  expect_equal(unname(m@bread), ref_bread, tolerance = 1e-5)

  # Check confidence intervals
  ci <- confidence_intervals(m)
  expect_equal(unname(ci[, "lower"]), ref$ci_lower, tolerance = 1e-5)
  expect_equal(unname(ci[, "upper"]), ref$ci_upper, tolerance = 1e-5)
})

test_that("ee_robust_regression Huber defaults to 'huber' loss", {
  ref <- load_fixture("ee_robust_regression_huber")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  k <- ref$k

  # The default loss should be "huber", so omitting it should give the same
  # results as explicitly specifying loss = "huber"
  psi <- function(theta) {
    ee_robust_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      k = k
    )
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m, solver = "lm")

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-5)
})

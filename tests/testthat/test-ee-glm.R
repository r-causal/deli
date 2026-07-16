# Tests for ee_glm() estimating equations using Python fixtures
# ee_glm() is a general GLM estimating equation that takes theta, X, y,
# distribution, and link arguments. For canonical links (identity/normal,
# logit/binomial, log/poisson), the deriv/variance ratio simplifies to 1,
# so the results should match ee_regression.

test_that("ee_glm normal/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with normal distribution and identity link is equivalent to linear
  # regression. The estimating equation simplifies to:
  # psi(theta) = (y - X %*% theta) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "normal", link = "identity")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm gaussian/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # "gaussian" is an alias for "normal" in the Python implementation
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gaussian", link = "identity")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm binomial/logit matches ee_regression logistic", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with binomial distribution and logit link is equivalent to logistic

  # regression. For this canonical link, deriv/variance = 1, so the
  # estimating equation simplifies to: psi(theta) = (y - expit(X %*% theta)) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "logit")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm binomial/logistic link alias matches ee_regression logistic", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # "logistic" is an alias for "logit" link in the Python implementation
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "logistic")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

test_that("ee_glm poisson/log matches ee_regression poisson", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # GLM with Poisson distribution and log link is the canonical Poisson
  # regression. For this canonical link, deriv/variance = 1, so the
  # estimating equation simplifies to: psi(theta) = (y - exp(X %*% theta)) * X
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "poisson", link = "log")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-3)
  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-3)
})

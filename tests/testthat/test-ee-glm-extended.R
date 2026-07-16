# Extended tests for ee_glm() estimating equations
# These tests cover non-canonical link/distribution combinations and verify
# equivalence with ee_regression for canonical links using fixture data.

test_that("ee_glm probit/binomial converges", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Probit link with binomial distribution is a non-canonical GLM.
  # We verify that estimation converges and produces reasonable results.
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "probit")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Should converge: theta should be finite

  expect_true(all(is.finite(m@theta)))
  # Variance should be positive on diagonal
  expect_true(all(diag(m@variance) > 0))
  # Theta length should match init
  expect_length(m@theta, length(init))
})

test_that("ee_glm cloglog/binomial converges", {
  ref <- load_fixture("ee_regression_logistic")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Complementary log-log link with binomial distribution is a non-canonical
  # GLM. We verify that estimation converges and produces reasonable results.
  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "binomial", link = "cloglog")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Should converge: theta should be finite
  expect_true(all(is.finite(m@theta)))
  # Variance should be positive on diagonal
  expect_true(all(diag(m@variance) > 0))
  # Theta length should match init
  expect_length(m@theta, length(init))
})

test_that("ee_glm poisson/log matches ee_regression poisson", {
  ref <- load_fixture("ee_regression_poisson")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Poisson/log is a canonical link, so ee_glm should give the same results
  # as ee_regression with model = "poisson"
  psi_glm <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "poisson", link = "log")
  }

  psi_reg <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "poisson")
  }

  m_glm <- MEstimator(stacked_equations = psi_glm, init = init)
  m_glm <- estimate(m_glm)

  m_reg <- MEstimator(stacked_equations = psi_reg, init = init)
  m_reg <- estimate(m_reg)

  expect_equal(m_glm@theta, m_reg@theta, tolerance = 1e-6)
  expect_equal(diag(m_glm@variance), diag(m_reg@variance), tolerance = 1e-6)
})

test_that("ee_glm normal/identity matches ee_regression linear", {
  ref <- load_fixture("ee_regression_linear")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Normal/identity is a canonical link, so ee_glm should give the same
  # results as ee_regression with model = "linear"
  psi_glm <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "normal", link = "identity")
  }

  psi_reg <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }

  m_glm <- MEstimator(stacked_equations = psi_glm, init = init)
  m_glm <- estimate(m_glm)

  m_reg <- MEstimator(stacked_equations = psi_reg, init = init)
  m_reg <- estimate(m_reg)

  expect_equal(m_glm@theta, m_reg@theta, tolerance = 1e-6)
  expect_equal(diag(m_glm@variance), diag(m_reg@variance), tolerance = 1e-6)
})

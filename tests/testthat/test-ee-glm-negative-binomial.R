# Tests for ee_glm() with the negative binomial distribution (bd-3cdf.1)
# The negative binomial GLM estimates the regression coefficients plus one
# additional dispersion parameter (on the log scale), so theta has length
# ncol(X) + 1. Python appends a polygamma-based nuisance estimating equation
# for that dispersion parameter and carries it honestly through the sandwich
# variance. These tests encode that behavior and are backed by fixtures
# generated from Python Delicatessen via
# generate_glm_negative_binomial_fixtures.py.

test_that("ee_glm negative_binomial appends the polygamma nuisance equation", {
  ref <- load_fixture("ee_glm_negative_binomial_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # The raw estimating function must return ncol(X) + 1 rows: one score row
  # per regression coefficient plus the dispersion nuisance row.
  ee <- ee_glm(init, X = X, y = y,
               distribution = "negative_binomial", link = "log")

  expect_equal(nrow(ee), ref$ee_nrow)
  expect_equal(nrow(ee), ncol(X) + 1L)
  expect_equal(ncol(ee), length(y))

  # Row sums at the starting values pin the exact form of every equation,
  # including the polygamma nuisance row for the dispersion parameter.
  expect_equal(rowSums(ee), ref$ee_sum_at_init, tolerance = 1e-6)
})

test_that("ee_glm negative_binomial/log matches Python Delicatessen", {
  ref <- load_fixture("ee_glm_negative_binomial_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y,
           distribution = "negative_binomial", link = "log")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # theta has length ncol(X) + 1 (regression coefficients plus log-dispersion).
  expect_length(m@theta, ncol(X) + 1L)

  # theta, meat, and bread match Python to ~1e-8, but the sandwich variance is
  # compared at 1e-5. The dispersion parameter's bread entry comes from
  # numerically differentiating the polygamma nuisance equation at dx = 1e-9,
  # which sits near floating-point precision; inverting the bread amplifies
  # that noise into the variance. This matches the tolerance the causal
  # sandwich fixtures already use for the same reason.
  expect_python_match(m, "ee_glm_negative_binomial_log", tolerance = 1e-5)
})

test_that("ee_glm 'nb' alias matches the full negative_binomial name", {
  ref <- load_fixture("ee_glm_nb_alias_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "nb", link = "log")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_length(m@theta, ncol(X) + 1L)
  # See the note above on the 1e-5 sandwich-variance tolerance for the
  # numerically differentiated dispersion nuisance equation.
  expect_python_match(m, "ee_glm_nb_alias_log", tolerance = 1e-5)
})

test_that("ee_glm negative_binomial leaves the nuisance row unweighted", {
  ref <- load_fixture("ee_glm_negative_binomial_weighted_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  weights <- ref$weights

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y,
           distribution = "negative_binomial", link = "log",
           weights = weights)
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # Python weights the beta score rows but deliberately leaves the dispersion
  # nuisance row unweighted (regression.py:341). This weighted fit reproduces
  # that reference; weighting the nuisance row would move the dispersion
  # estimate away from Python. See the note above on the 1e-5 tolerance for the
  # numerically differentiated nuisance equation.
  expect_length(m@theta, ncol(X) + 1L)
  expect_python_match(m, "ee_glm_negative_binomial_weighted_log",
                      tolerance = 1e-5)
})

# Tests for ee_glm() with the gamma distribution (bd-3cdf.4)
# The gamma GLM estimates the regression coefficients plus one additional shape
# parameter alpha (on the log scale), so theta has length ncol(X) + 1. The
# reciprocal of alpha is the McCullagh and Nelder GLM dispersion (phi = 1 /
# alpha). Python appends a digamma-based nuisance estimating equation for that
# shape parameter and carries it honestly through the sandwich variance. These
# tests encode that behavior and are backed by fixtures generated from Python
# Delicatessen via generate_glm_gamma_fixtures.py.

test_that("ee_glm gamma appends the shape nuisance equation", {
  ref <- load_fixture("ee_glm_gamma_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # The raw estimating function must return ncol(X) + 1 rows: one score row per
  # regression coefficient plus the shape nuisance row.
  ee <- ee_glm(init, X = X, y = y, distribution = "gamma", link = "log")

  expect_equal(nrow(ee), ref$ee_nrow)
  expect_equal(nrow(ee), ncol(X) + 1L)
  expect_equal(ncol(ee), length(y))

  # Row sums at the starting values pin the exact form of every equation,
  # including the shape nuisance row.
  expect_equal(rowSums(ee), ref$ee_sum_at_init, tolerance = 1e-6)
})

test_that("ee_glm gamma/log matches Python Delicatessen", {
  ref <- load_fixture("ee_glm_gamma_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gamma", link = "log")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # theta has length ncol(X) + 1 (regression coefficients plus log-shape).
  expect_length(m@theta, ncol(X) + 1L)

  # theta, bread, meat, and the full sandwich variance all reproduce Python
  # within the default 1e-6 tolerance.
  expect_python_match(m, "ee_glm_gamma_log", tolerance = 1e-6)
})

test_that("ee_glm gamma/inverse (non-canonical link) matches Python", {
  ref <- load_fixture("ee_glm_gamma_inverse")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gamma", link = "inverse")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_length(m@theta, ncol(X) + 1L)
  expect_python_match(m, "ee_glm_gamma_inverse", tolerance = 1e-6)
})

test_that("ee_glm gamma weights the shape nuisance row like Python", {
  ref <- load_fixture("ee_glm_gamma_weighted_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  weights <- ref$weights

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gamma", link = "log",
           weights = weights)
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # The weighted fit must reproduce Python, whose gamma nuisance row is weighted.
  # Dropping the weights from the nuisance row would shift the shape estimate
  # away from this reference.
  expect_length(m@theta, ncol(X) + 1L)
  expect_python_match(m, "ee_glm_gamma_weighted_log", tolerance = 1e-6)

  # Independent cross-check: because the weights are integers, the weighted fit
  # must equal the unweighted fit on the row-expanded (repeated) data. This
  # holds for all ncol(X) + 1 parameters only if the shape nuisance row carries
  # the weights; an unweighted nuisance would break the shape equality.
  idx <- rep(seq_len(nrow(X)), times = weights)
  X_expanded <- X[idx, , drop = FALSE]
  y_expanded <- y[idx]

  psi_expanded <- function(theta) {
    ee_glm(theta, X = X_expanded, y = y_expanded,
           distribution = "gamma", link = "log")
  }

  m_expanded <- estimate(MEstimator(stacked_equations = psi_expanded,
                                    init = init))

  expect_equal(unname(m@theta), unname(m_expanded@theta), tolerance = 1e-6)
})

# Tests for ee_glm() with the gamma distribution
# The gamma GLM estimates the regression coefficients plus one additional shape
# parameter alpha (on the log scale), so theta has length ncol(X) + 1. The
# reciprocal of alpha is the McCullagh and Nelder GLM dispersion (phi = 1 /
# alpha). Python appends a digamma-based nuisance estimating equation for that
# shape parameter and carries it honestly through the sandwich variance. These
# tests encode that behavior and are backed by fixtures generated from Python
# Delicatessen via generate_glm_gamma_fixtures.py.
#
# Tolerance posture: the shape nuisance row is differentiated numerically, and
# the sandwich amplifies that derivative's error, so the variance and asymptotic
# variance blocks move between platforms in a way the rest of the comparison
# does not. The fixture comparisons below keep theta, the bread, and the
# interval bounds pinned at 1e-6 and give the variance blocks 1e-4. Minimum
# passing tolerances measured on macOS and on two Linux images (R 4.4.3 and
# R 4.6): theta at or under 2.5e-11, bread at or under 1.3e-7, interval bounds
# at or under 1.1e-7, and the variance blocks between 1.2e-7 and 8.7e-7. CI has
# seen the variance blocks pass 1e-6 on Windows and on Ubuntu oldrel-2, with
# element-wise relative differences reaching 4e-5 on the small off-diagonal
# entries, which 1e-4 covers with room left.

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
  # including the shape nuisance row. The Python reference carries no labels,
  # so the comparison is made on the values alone, as the fixture helper does.
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)
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

  # theta, bread, and the interval bounds reproduce Python at 1e-6. The variance
  # blocks measure 2.9e-7 on macOS and 7.8e-7 on both Linux images; this is the
  # comparison the Ubuntu oldrel-2 CI job pushed past 1e-6.
  expect_python_match(
    m,
    "ee_glm_gamma_log",
    tolerance = 1e-6,
    variance_tolerance = 1e-4
  )
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
  # The variance blocks measure 8.6e-7 on macOS, the widest of the gamma
  # fixtures, against 1.2e-7 on both Linux images; this is the comparison the
  # Windows CI job pushed past 1e-6.
  expect_python_match(
    m,
    "ee_glm_gamma_inverse",
    tolerance = 1e-6,
    variance_tolerance = 1e-4
  )
})

test_that("ee_glm gamma weights the shape nuisance row like Python", {
  ref <- load_fixture("ee_glm_gamma_weighted_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init
  weights <- ref$weights

  psi <- function(theta) {
    ee_glm(
      theta,
      X = X,
      y = y,
      distribution = "gamma",
      link = "log",
      weights = weights
    )
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # The weighted fit must reproduce Python, whose gamma nuisance row is weighted.
  # Dropping the weights from the nuisance row would shift the shape estimate
  # away from this reference.
  expect_length(m@theta, ncol(X) + 1L)
  # The variance blocks measure 7.8e-7 on macOS, 3.1e-7 on R 4.4.3 Linux, and
  # 5.1e-7 on R 4.6 Linux. CI has not failed this one, but its margin sits with
  # the other two gamma fixtures, so it gets the same room.
  expect_python_match(
    m,
    "ee_glm_gamma_weighted_log",
    tolerance = 1e-6,
    variance_tolerance = 1e-4
  )

  # Independent cross-check: because the weights are integers, the weighted fit
  # must equal the unweighted fit on the row-expanded (repeated) data. This
  # holds for all ncol(X) + 1 parameters only if the shape nuisance row carries
  # the weights; an unweighted nuisance would break the shape equality.
  idx <- rep(seq_len(nrow(X)), times = weights)
  X_expanded <- X[idx, , drop = FALSE]
  y_expanded <- y[idx]

  psi_expanded <- function(theta) {
    ee_glm(
      theta,
      X = X_expanded,
      y = y_expanded,
      distribution = "gamma",
      link = "log"
    )
  }

  m_expanded <- estimate(MEstimator(
    stacked_equations = psi_expanded,
    init = init
  ))

  expect_equal(unname(m@theta), unname(m_expanded@theta), tolerance = 1e-6)
})

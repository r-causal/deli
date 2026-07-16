# Tests for ee_glm() with the tweedie distribution (bd-3cdf.3)
#
# The tweedie distribution generalizes the Poisson and gamma GLMs through a
# fixed (not estimated) power hyperparameter p that sets the variance function
# v(mu) = mu^p. Unlike the negative binomial and gamma distributions, tweedie
# does NOT estimate an additional parameter, so theta has length ncol(X) and no
# nuisance estimating equation is appended. Notable equivalences: p = 1 is the
# Poisson GLM, and p = 2 shares the gamma GLM's beta score equations.
#
# These tests are backed by fixtures generated from Python Delicatessen via
# generate_glm_tweedie_fixtures.py. Comparisons use the default 1e-6 tolerance:
# tweedie estimating equations are pure beta scores (no polygamma nuisance row
# like the negative binomial), so the numerically differentiated bread is smooth
# and the sandwich variance reproduces Python well within 1e-6.

test_that("ee_glm tweedie returns ncol(X) rows with no nuisance equation", {
  ref <- load_fixture("ee_glm_tweedie_p15_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # The raw estimating function must return exactly ncol(X) rows: one score row
  # per regression coefficient and nothing more. This discriminates tweedie from
  # the negative binomial and gamma distributions, which append a nuisance row.
  ee <- ee_glm(init, X = X, y = y,
               distribution = "tweedie", link = "log", hyperparameter = 1.5)

  expect_equal(nrow(ee), ref$ee_nrow)
  expect_equal(nrow(ee), ncol(X))
  expect_equal(ncol(ee), length(y))

  # Row sums at the starting values pin the exact form of every score equation.
  expect_equal(rowSums(ee), ref$ee_sum_at_init, tolerance = 1e-6)
})

test_that("ee_glm tweedie hyperparameter=1.5 matches Python Delicatessen", {
  ref <- load_fixture("ee_glm_tweedie_p15_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_glm(theta, X = X, y = y,
           distribution = "tweedie", link = "log", hyperparameter = 1.5)
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # theta has length ncol(X); tweedie estimates no extra parameter.
  expect_length(m@theta, ncol(X))
  expect_python_match(m, "ee_glm_tweedie_p15_log", tolerance = 1e-6)
})

test_that("ee_glm tweedie hyperparameter=1 reproduces the Poisson GLM", {
  ref <- load_fixture("ee_glm_tweedie_p1_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  tweedie_psi <- function(theta) {
    ee_glm(theta, X = X, y = y,
           distribution = "tweedie", link = "log", hyperparameter = 1.0)
  }
  m_tw <- estimate(MEstimator(stacked_equations = tweedie_psi, init = init))

  expect_length(m_tw@theta, ncol(X))
  expect_python_match(m_tw, "ee_glm_tweedie_p1_log", tolerance = 1e-6)

  # At p = 1 the tweedie variance mu^1 equals the Poisson variance mu, so the
  # tweedie fit must coincide with the Poisson GLM on the same data (theta and
  # the full sandwich variance).
  poisson_psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "poisson", link = "log")
  }
  m_pois <- estimate(MEstimator(stacked_equations = poisson_psi, init = init))

  expect_equal(unname(m_tw@theta), unname(m_pois@theta), tolerance = 1e-6)
  expect_equal(unname(m_tw@variance), unname(m_pois@variance), tolerance = 1e-6)
})

test_that("ee_glm tweedie hyperparameter=2 matches Python and the gamma GLM beta", {
  ref <- load_fixture("ee_glm_tweedie_p2_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  tweedie_psi <- function(theta) {
    ee_glm(theta, X = X, y = y,
           distribution = "tweedie", link = "log", hyperparameter = 2.0)
  }
  m_tw <- estimate(MEstimator(stacked_equations = tweedie_psi, init = init))

  expect_length(m_tw@theta, ncol(X))
  expect_python_match(m_tw, "ee_glm_tweedie_p2_log", tolerance = 1e-6)

  # At p = 2 the tweedie variance mu^2 equals the gamma variance, so the beta
  # score equations coincide. The gamma GLM additionally estimates a shape
  # nuisance (its theta has one extra element), but its beta estimates must
  # match the tweedie p = 2 estimates on the same strictly positive data.
  gamma_psi <- function(theta) {
    ee_glm(theta, X = X, y = y, distribution = "gamma", link = "log")
  }
  m_gamma <- estimate(MEstimator(stacked_equations = gamma_psi,
                                 init = c(init, 0)))

  expect_length(m_gamma@theta, ncol(X) + 1L)
  expect_equal(unname(m_tw@theta),
               unname(m_gamma@theta)[seq_len(ncol(X))],
               tolerance = 1e-6)
})

test_that("ee_glm tweedie hyperparameter changes the solution", {
  # Same data, two different powers. If the hyperparameter were ignored the two
  # fits would be identical; they must differ, and each must match its fixture.
  ref1 <- load_fixture("ee_glm_tweedie_p1_log")
  ref15 <- load_fixture("ee_glm_tweedie_p15_log")

  X <- ref1$X
  y <- ref1$y
  init <- ref1$init

  fit <- function(p) {
    estimate(MEstimator(
      stacked_equations = function(theta) {
        ee_glm(theta, X = X, y = y,
               distribution = "tweedie", link = "log", hyperparameter = p)
      },
      init = init
    ))
  }

  m_p1 <- fit(1.0)
  m_p15 <- fit(1.5)

  # The two solutions are genuinely different (not a tolerance-floor artifact):
  # the largest coefficient difference is well above any solver noise.
  expect_gt(max(abs(unname(m_p1@theta) - unname(m_p15@theta))), 1e-3)

  # And each fit lands on its own Python reference.
  expect_equal(unname(m_p1@theta), ref1$theta, tolerance = 1e-6)
  expect_equal(unname(m_p15@theta), ref15$theta, tolerance = 1e-6)
})

test_that("ee_glm tweedie rejects a negative hyperparameter", {
  ref <- load_fixture("ee_glm_tweedie_p15_log")

  X <- ref$X
  y <- ref$y
  init <- ref$init

  # Python raises when the tweedie hyperparameter is negative (it must be
  # non-negative). The R implementation must mirror this, and the message must
  # name the "non-negative" constraint so this assertion does not pass on the
  # generic "distribution not supported" or "unused argument" errors.
  expect_error(
    ee_glm(init, X = X, y = y,
           distribution = "tweedie", link = "log", hyperparameter = -1),
    "non-negative"
  )
})

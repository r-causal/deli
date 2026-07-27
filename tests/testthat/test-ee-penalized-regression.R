# Tests for penalized regression EEs (bridge, LASSO, dLASSO, elastic net)

test_that("ee_bridge_regression with gamma=2 matches ee_ridge_regression", {
  ref <- load_fixture("ee_ridge_regression")
  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty

  psi_bridge <- function(theta) {
    ee_bridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      gamma = 2
    )
  }
  psi_ridge <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty
    )
  }

  m_bridge <- estimate(MEstimator(stacked_equations = psi_bridge, init = init))
  m_ridge <- estimate(MEstimator(stacked_equations = psi_ridge, init = init))

  expect_equal(m_bridge@theta, m_ridge@theta, tolerance = 1e-6)
  expect_equal(
    diag(m_bridge@variance),
    diag(m_ridge@variance),
    tolerance = 1e-6
  )
})

# The LASSO and elastic net reach bridge_penalty with gamma just above 1, so
# every evaluation of their estimating equations warns that the penalty is not
# everywhere differentiable and the sandwich variance is therefore not well
# defined. That warning is asserted on its own below; the tests that solve or
# compare fits suppress it at the call site so the assertion under test is not
# buried under repeated copies of the same warning, one per EE evaluation.

test_that("ee_lasso_regression returns correct shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  result <- suppressWarnings(
    ee_lasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 0.5
    )
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(X))
  expect_equal(ncol(result), nrow(X))
})

test_that("ee_lasso_regression solves via MEstimator", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    suppressWarnings(
      ee_lasso_regression(theta, X = X, y = y, model = "linear", penalty = 0.5)
    )
  }
  m <- suppressWarnings(estimate(MEstimator(
    stacked_equations = psi,
    init = init
  )))
  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_lasso_regression with zero penalty matches unpenalized", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  # The warning depends only on gamma, so a zero penalty still raises it
  psi_lasso <- function(theta) {
    suppressWarnings(
      ee_lasso_regression(theta, X = X, y = y, model = "linear", penalty = 0)
    )
  }
  psi_linear <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }

  m_lasso <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi_lasso, init = init))
  )
  m_linear <- estimate(MEstimator(stacked_equations = psi_linear, init = init))

  expect_equal(m_lasso@theta, m_linear@theta, tolerance = 1e-4)
})

test_that("ee_dlasso_regression returns correct shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  result <- ee_dlasso_regression(
    init,
    X = X,
    y = y,
    model = "linear",
    penalty = 0.5
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(X))
})

test_that("ee_dlasso_regression solves via MEstimator", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi <- function(theta) {
    ee_dlasso_regression(theta, X = X, y = y, model = "linear", penalty = 0.5)
  }
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  expect_true(all(is.finite(m@theta)))
})

test_that("ee_elasticnet_regression returns correct shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  result <- suppressWarnings(
    ee_elasticnet_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 0.5,
      ratio = 0.5
    )
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), ncol(X))
})

test_that("ee_elasticnet_regression with ratio=0 matches ridge", {
  ref <- load_fixture("ee_ridge_regression")
  X <- ref$X
  y <- ref$y
  init <- ref$init
  penalty <- ref$penalty

  # ratio = 0 zeroes the L1 contribution, but the L1 penalty is still evaluated
  # to form it, so the non-differentiability warning is still raised
  psi_en <- function(theta) {
    suppressWarnings(
      ee_elasticnet_regression(
        theta,
        X = X,
        y = y,
        model = "linear",
        penalty = penalty,
        ratio = 0
      )
    )
  }
  psi_ridge <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty
    )
  }

  m_en <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi_en, init = init))
  )
  m_ridge <- estimate(MEstimator(stacked_equations = psi_ridge, init = init))

  expect_equal(m_en@theta, m_ridge@theta, tolerance = 1e-4)
})

test_that("ee_elasticnet_regression errors for invalid ratio", {
  ref <- load_fixture("ee_regression_linear")
  expect_error(
    ee_elasticnet_regression(
      ref$init,
      X = ref$X,
      y = ref$y,
      model = "linear",
      penalty = 0.5,
      ratio = 1.5
    ),
    "ratio"
  )
})

# Behaviour-preservation pins for the performance refactors -------------------
#
# The penalized estimating functions weight the whole equation, penalty
# included, and recycle a length-p penalty down the rows of a p-by-n matrix.
# These reconstruct that matrix element by element from the definition, so a
# refactor of the transpose and weighting structure is caught to machine
# precision, including the non-unit-weight and offset paths that make the
# weighting placement observable.

manual_penalized <- function(theta, X, y, penalty_terms, weights, offset) {
  n <- nrow(X)
  p <- ncol(X)
  eta <- as.numeric(X %*% theta) + offset
  resid <- y - eta # linear model
  out <- matrix(0, nrow = p, ncol = n)
  for (i in seq_len(n)) {
    out[, i] <- (X[i, ] * resid[i] - penalty_terms) * weights[i]
  }
  out
}

bridge_penalty_terms <- function(theta, n, penalty, gamma, center = 0) {
  penalty <- rep(as.numeric(penalty), length.out = length(theta))
  center <- rep(as.numeric(center), length.out = length(theta))
  (penalty / (gamma * n)) *
    gamma *
    abs(theta - center)^(gamma - 1) *
    sign(theta - center)
}

test_that("ee_bridge_regression matrix matches the element-wise definition with weights and offset", {
  set.seed(21)
  n <- 50
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- as.numeric(X %*% c(1, -0.5, 0.3)) + rnorm(n)
  theta <- c(0.4, 0.6, -0.2)
  penalty <- c(0, 8, 8)
  gamma <- 2
  weights <- runif(n, 0.3, 2.5)
  offset <- rnorm(n, sd = 0.3)

  pt <- bridge_penalty_terms(theta, n, penalty, gamma)
  expect_equal(
    ee_bridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      gamma = gamma,
      weights = weights,
      offset = offset
    ),
    manual_penalized(theta, X, y, pt, weights, offset),
    tolerance = 1e-12
  )
})

test_that("ee_dlasso_regression and ee_elasticnet_regression weight the penalty inside the transpose", {
  set.seed(22)
  n <- 45
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- as.numeric(X %*% c(0.5, 1, -0.4)) + rnorm(n)
  theta <- c(0.3, 0.7, -0.1)
  penalty <- c(0, 6, 6)
  weights <- runif(n, 0.4, 2)
  offset <- rnorm(n, sd = 0.25)

  # dLASSO penalty terms (standard-normal smoothing of the L1 penalty)
  tc <- theta - 0
  s <- 1e-6
  pt_dlasso <- (penalty / n) *
    (2 * pnorm(tc / s) + 2 * (tc / s) * dnorm(tc / s) - 1)
  expect_equal(
    ee_dlasso_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = weights,
      offset = offset
    ),
    manual_penalized(theta, X, y, pt_dlasso, weights, offset),
    tolerance = 1e-10
  )

  # elastic net combines an approximate-LASSO bridge (gamma = 1 + epsilon) with a
  # ridge bridge (gamma = 2) at the given ratio
  ratio <- 0.5
  epsilon <- 3e-3
  pt_en <- ratio *
    bridge_penalty_terms(theta, n, penalty, 1 + epsilon) +
    (1 - ratio) * bridge_penalty_terms(theta, n, penalty, 2)
  expect_equal(
    suppressWarnings(ee_elasticnet_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      ratio = ratio,
      weights = weights,
      offset = offset
    )),
    manual_penalized(theta, X, y, pt_en, weights, offset),
    tolerance = 1e-10
  )
})

test_that("NULL weights and unit weights give identical psi for ridge and dlasso", {
  set.seed(24)
  n <- 40
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- as.numeric(X %*% c(0.7, -0.3, 0.5)) + rnorm(n)
  theta <- c(0.2, 0.4, -0.1)
  penalty <- c(0, 7, 7)

  expect_identical(
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = NULL
    ),
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = rep(1, n)
    )
  )

  expect_identical(
    ee_dlasso_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = NULL
    ),
    ee_dlasso_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = rep(1, n)
    )
  )
})

test_that("ee_ridge_regression weighted exact bread agrees with the central difference", {
  set.seed(23)
  n <- 120
  X <- cbind(1, rnorm(n), rnorm(n))
  y <- as.numeric(X %*% c(0.8, -0.6, 0.4)) + rnorm(n)
  weights <- runif(n, 0.5, 2)
  penalty <- c(0, 5, 5)

  psi <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty,
      weights = weights
    )
  }
  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = c(0, 0, 0)),
    deriv_method = "exact"
  )
  m_cap <- estimate(
    MEstimator(stacked_equations = psi, init = c(0, 0, 0)),
    deriv_method = "capprox"
  )
  expect_equal(unname(m_exact@bread), unname(m_cap@bread), tolerance = 1e-6)
  expect_equal(
    sqrt(diag(m_exact@variance)),
    sqrt(diag(m_cap@variance)),
    tolerance = 1e-6
  )
})

# Input validation (mirrors Python check_penalty_shape and the bridge/dlasso
# hyperparameter checks). These pin that validation.
#
# Python messages to mirror (delicatessen/errors.py, regression.py):
#   penalty shape : "The penalty term must be either a single number or the
#                    same length as theta."
#   center shape  : "The center term must be either a single number or the
#                    same length as theta."
#   negative      : "All penalty terms must be non-negative"
#   gamma < 1     : "L_{gamma} for `gamma` < 1 is not currently able to be
#                    supported with estimating equations evaluated using
#                    numerical methods."
#   gamma < 2     : (UserWarning) "The estimating equation for chosen penalized
#                    regression model is not always differentiable. Therefore,
#                    the bread matrix is not always defined for finite samples,
#                    and the sandwich should not be used to estimate the
#                    variance."
#   dlasso s < 0  : "`s` must be greater than zero for the approximate LASSO"
#
# Each call is valid except for the one input under test, so the
# implementation raises only because of that input. The EEs are called
# directly (not through the solver) so the raised condition comes from the
# validation, not from an incidental solver error.

test_that("ee_ridge_regression validates penalty and center shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init # length 2

  # Penalty neither scalar nor length(theta)
  expect_error(
    ee_ridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(1, 1, 1)
    ),
    "single number or the same"
  )

  # Negative penalty
  expect_error(
    ee_ridge_regression(init, X = X, y = y, model = "linear", penalty = -1),
    "non-negative"
  )

  # Center neither scalar nor length(theta); penalty kept valid
  expect_error(
    ee_ridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      center = c(0, 0, 0)
    ),
    "center"
  )
})

test_that("ee_bridge_regression validates penalty and center shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(1, 1, 1),
      gamma = 2
    ),
    "single number or the same"
  )

  # Element-wise non-negativity: one entry negative
  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(0.5, -1),
      gamma = 2
    ),
    "non-negative"
  )

  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      gamma = 2,
      center = c(0, 0, 0)
    ),
    "center"
  )
})

test_that("ee_bridge_regression errors for gamma below 1", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      gamma = 0.5
    ),
    "gamma"
  )
})

test_that("ee_bridge_regression warns for gamma in [1, 2)", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_warning(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      gamma = 1.5
    ),
    "differentiable"
  )
})

test_that("ee_bridge_regression does not warn for gamma >= 2", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  # gamma = 2 (ridge) is differentiable; no non-differentiability warning
  expect_no_warning(
    ee_bridge_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      gamma = 2
    )
  )
})

test_that("ee_lasso_regression warns that it is not always differentiable", {
  ref <- load_fixture("ee_regression_linear")

  expect_warning(
    ee_lasso_regression(
      ref$init,
      X = ref$X,
      y = ref$y,
      model = "linear",
      penalty = 1
    ),
    "not always differentiable"
  )
})

test_that("ee_elasticnet_regression warns that it is not always differentiable", {
  ref <- load_fixture("ee_regression_linear")

  expect_warning(
    ee_elasticnet_regression(
      ref$init,
      X = ref$X,
      y = ref$y,
      model = "linear",
      penalty = 1,
      ratio = 0.5
    ),
    "not always differentiable"
  )
})

test_that("ee_lasso_regression validates penalty shape and sign", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_lasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(1, 1, 1)
    ),
    "single number or the same"
  )

  expect_error(
    ee_lasso_regression(init, X = X, y = y, model = "linear", penalty = -1),
    "non-negative"
  )
})

test_that("ee_dlasso_regression validates penalty and center shape", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_dlasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(1, 1, 1)
    ),
    "single number or the same"
  )

  expect_error(
    ee_dlasso_regression(init, X = X, y = y, model = "linear", penalty = -1),
    "non-negative"
  )

  expect_error(
    ee_dlasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      center = c(0, 0, 0)
    ),
    "center"
  )
})

test_that("ee_dlasso_regression errors for s below zero", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_dlasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      s = -1
    ),
    "greater than zero"
  )
})

test_that("ee_elasticnet_regression validates penalty shape and sign", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_elasticnet_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = c(1, 1, 1),
      ratio = 0.5
    ),
    "single number or the same"
  )

  expect_error(
    ee_elasticnet_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = -1,
      ratio = 0.5
    ),
    "non-negative"
  )
})

test_that("stronger penalty shrinks coefficients more", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  psi_low <- function(theta) {
    suppressWarnings(
      ee_lasso_regression(theta, X = X, y = y, model = "linear", penalty = 0.1)
    )
  }
  psi_high <- function(theta) {
    suppressWarnings(
      ee_lasso_regression(theta, X = X, y = y, model = "linear", penalty = 10)
    )
  }

  m_low <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi_low, init = init))
  )
  m_high <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi_high, init = init))
  )

  # Higher penalty should shrink coefficients more (smaller L2 norm)
  expect_true(
    sum(m_high@theta^2) <= sum(m_low@theta^2) + 1e-6
  )
})

# Weighted-penalty parity.
#
# Python weights the entire penalized estimating equation, penalty included:
#
#   w * (((y - pred_y) * X).T - penalty_terms[:, None])
#
# An implementation that subtracts the penalty outside the weight differs from
# this by penalty_term * (sum(weights) - n) in each estimating-function row sum,
# which propagates into theta, variance, and bread once the equation is solved.
# With unit weights sum(weights) equals n and the placement is invisible, so
# non-unit weights are what make these tests bite.
#
# Fixtures: generate_weighted_penalized_fixtures.py (seed 20240706, n = 100,
# integer weights in {1,...,5} with sum(weights) = 308, penalty = c(0, 10, 10)
# leaving the intercept unpenalized). The init c(0.5, 0.5, -0.3) is deliberately
# non-zero so the penalty term is non-zero at init and the raw-EE pin below
# discriminates without a solve.
#
# On these fixtures the unweighted-penalty form is orders of magnitude away from
# Python across every model: theta gaps of 0.047 to 0.18 (largest on the
# additive spline terms) and EE row-sum gaps of 2.9 to 20.8, all far above the
# 1e-5 parity and 1e-6 raw-EE tolerances used below.
#
# Parity tolerance is 1e-5, matching the unweighted penalized fixture tests in
# test-python-fixtures.R (the penalized solves do not settle to 1e-6). The
# raw-EE-at-init pin uses 1e-6 because it is a direct evaluation, no solve.
#
# The approximate LASSO and elastic net are not everywhere differentiable and
# warn on every evaluation; those call sites are wrapped in suppressWarnings so
# the parity assertions stay uncluttered.

test_that("ee_ridge_regression weighted matches Python (penalty inside weight)", {
  ref <- load_fixture("ee_ridge_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  # Raw estimating-function row sums at init pin the weighted-penalty form
  # directly, without a solve.
  ee <- ee_ridge_regression(
    init,
    X = X,
    y = y,
    model = "linear",
    penalty = ref$penalty,
    weights = w
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      weights = w
    )
  }
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  expect_python_match(m, "ee_ridge_regression_weighted", tolerance = 1e-5)
})

test_that("ee_bridge_regression weighted matches Python (penalty inside weight)", {
  ref <- load_fixture("ee_bridge_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  ee <- ee_bridge_regression(
    init,
    X = X,
    y = y,
    model = "linear",
    penalty = ref$penalty,
    gamma = ref$gamma,
    weights = w
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_bridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      gamma = ref$gamma,
      weights = w
    )
  }
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  expect_python_match(m, "ee_bridge_regression_weighted", tolerance = 1e-5)
})

test_that("ee_lasso_regression weighted matches Python (penalty inside weight)", {
  ref <- load_fixture("ee_lasso_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  ee <- suppressWarnings(
    ee_lasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      weights = w
    )
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_lasso_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      weights = w
    )
  }
  m <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi, init = init))
  )
  expect_python_match(m, "ee_lasso_regression_weighted", tolerance = 1e-5)
})

test_that("ee_dlasso_regression weighted matches Python (penalty inside weight)", {
  ref <- load_fixture("ee_dlasso_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  ee <- ee_dlasso_regression(
    init,
    X = X,
    y = y,
    model = "linear",
    penalty = ref$penalty,
    weights = w
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_dlasso_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      weights = w
    )
  }
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  expect_python_match(m, "ee_dlasso_regression_weighted", tolerance = 1e-5)
})

test_that("ee_elasticnet_regression weighted matches Python (penalty inside weight)", {
  ref <- load_fixture("ee_elasticnet_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  ee <- suppressWarnings(
    ee_elasticnet_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      ratio = ref$ratio,
      weights = w
    )
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_elasticnet_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = ref$penalty,
      ratio = ref$ratio,
      weights = w
    )
  }
  m <- suppressWarnings(
    estimate(MEstimator(stacked_equations = psi, init = init))
  )
  expect_python_match(m, "ee_elasticnet_regression_weighted", tolerance = 1e-5)
})

test_that("ee_additive_regression weighted matches Python (delegates to bridge)", {
  ref <- load_fixture("ee_additive_regression_weighted")
  X <- ref$X
  y <- ref$y
  w <- ref$weights
  init <- ref$init

  # Reconstruct specifications from the fixture; jsonlite mangles nested list
  # structures under simplifyVector = TRUE (see the unweighted additive test).
  specs <- list(
    NULL,
    list(
      knots = c(-2, -1, 0, 1, 2),
      penalty = 5,
      natural = TRUE,
      power = 3,
      normalized = FALSE
    )
  )

  ee <- ee_additive_regression(
    init,
    X = X,
    y = y,
    specifications = specs,
    model = "linear",
    weights = w
  )
  expect_equal(unname(rowSums(ee)), ref$ee_sum_at_init, tolerance = 1e-6)

  psi <- function(theta) {
    ee_additive_regression(
      theta,
      X = X,
      y = y,
      specifications = specs,
      model = "linear",
      weights = w
    )
  }
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  # The penalized-spline design yields an ill-conditioned bread whose numerical
  # derivative inversion noise differs across derivative implementations by more
  # than 1e-5, so the variance and asymptotic variance are compared at 1e-4,
  # matching the unweighted additive fixture test.
  expect_python_match(m, "ee_additive_regression_weighted", tolerance = 1e-4)
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_bridge_regression rejects a y whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = ref$y[seq_len(nrow(X) / 2)],
      model = "linear",
      penalty = 1,
      gamma = 2
    ),
    "same length as the data"
  )
})

test_that("ee_bridge_regression rejects an offset whose length differs from the rows of X", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  init <- ref$init

  expect_error(
    ee_bridge_regression(
      init,
      X = X,
      y = ref$y,
      model = "linear",
      penalty = 1,
      gamma = 2,
      offset = rep(0, nrow(X) / 2)
    ),
    "same length as the data"
  )
})

test_that("ee_lasso_regression rejects a negative epsilon with an epsilon-named message", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_lasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      epsilon = -0.5
    ),
    "epsilon"
  )
})

test_that("ee_elasticnet_regression rejects a negative epsilon with an epsilon-named message", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_elasticnet_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      ratio = 0.5,
      epsilon = -0.5
    ),
    "epsilon"
  )
})

test_that("ee_dlasso_regression rejects s = 0 to match its documentation", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y
  init <- ref$init

  expect_error(
    ee_dlasso_regression(
      init,
      X = X,
      y = y,
      model = "linear",
      penalty = 1,
      s = 0
    ),
    "greater than zero"
  )
})

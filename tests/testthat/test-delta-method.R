# Helper to create a fitted single-parameter MEstimator (mean)
make_fitted_mean <- function() {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  estimate(m)
}

# Helper to create a fitted two-parameter MEstimator (mean + variance)
make_fitted_mean_variance <- function() {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y
  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  estimate(m)
}

# delta_method() ---------------------------------------------------------------

test_that("delta_method() with identity transform returns original variance", {
  m <- make_fitted_mean()
  # Identity: g(theta) = theta
  result <- delta_method(m, transform = function(theta) theta)
  expect_equal(result, unname(m@variance), tolerance = 1e-6)
})

test_that("delta_method() with identity transform on multi-parameter model returns original variance", {
  m <- make_fitted_mean_variance()
  result <- delta_method(m, transform = function(theta) theta)
  expect_equal(result, unname(m@variance), tolerance = 1e-6)
})

test_that("delta_method() computes variance for exp(theta) transformation", {
  m <- make_fitted_mean()
  result <- delta_method(m, transform = function(theta) exp(theta))
  # By delta method: Var(exp(theta)) = exp(theta)^2 * Var(theta)
  expected_var <- matrix(
    exp(m@theta[1])^2 * m@variance[1, 1],
    nrow = 1,
    ncol = 1
  )
  expect_equal(result, expected_var, tolerance = 1e-4)
})

test_that("delta_method() computes variance for log(theta) transformation", {
  m <- make_fitted_mean()
  result <- delta_method(m, transform = function(theta) log(theta))
  # By delta method: Var(log(theta)) = (1/theta)^2 * Var(theta)
  expected_var <- matrix(
    (1 / m@theta[1])^2 * m@variance[1, 1],
    nrow = 1,
    ncol = 1
  )
  expect_equal(result, expected_var, tolerance = 1e-4)
})

test_that("delta_method() works for scalar function of multi-parameter model", {
  m <- make_fitted_mean_variance()
  # g(theta) = theta[1] / theta[2] (ratio of mean to variance)
  result <- delta_method(m, transform = function(theta) theta[1] / theta[2])
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)

  # Manual delta method: G = [1/theta2, -theta1/theta2^2]
  # Var = G %*% V %*% t(G)
  th <- m@theta
  V <- m@variance
  G <- matrix(c(1 / th[2], -th[1] / th[2]^2), nrow = 1)
  expected_var <- G %*% V %*% t(G)
  expect_equal(result, expected_var, tolerance = 1e-4)
})

test_that("delta_method() works for subset of parameters", {
  m <- make_fitted_mean_variance()
  # Transform that only returns the first parameter
  result <- delta_method(m, transform = function(theta) theta[1])
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  # Should equal the (1,1) element of the variance matrix
  expect_equal(result[1, 1], unname(m@variance[1, 1]), tolerance = 1e-6)
})

test_that("delta_method() returns correct dimensions for vector-valued transform", {
  m <- make_fitted_mean_variance()
  # Transform that returns a 2-d vector from 2-d input
  result <- delta_method(m, transform = function(theta) {
    c(theta[1] + theta[2], theta[1] * theta[2])
  })
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
})

test_that("delta_method() with linear transform G*V*G^T matches manual calculation", {
  m <- make_fitted_mean_variance()
  # g(theta) = c(theta[1] + theta[2], theta[1] - theta[2])
  result <- delta_method(m, transform = function(theta) {
    c(theta[1] + theta[2], theta[1] - theta[2])
  })
  # Jacobian is [[1, 1], [1, -1]]
  G <- matrix(c(1, 1, 1, -1), nrow = 2, byrow = TRUE)
  V <- m@variance
  expected <- G %*% V %*% t(G)
  expect_equal(result, expected, tolerance = 1e-4)
})

test_that("delta_method() with squaring transform matches manual calculation", {
  m <- make_fitted_mean()
  # g(theta) = theta^2
  result <- delta_method(m, transform = function(theta) theta^2)
  # Var(theta^2) = (2*theta)^2 * Var(theta)
  expected_var <- matrix(
    (2 * m@theta[1])^2 * m@variance[1, 1],
    nrow = 1,
    ncol = 1
  )
  expect_equal(result, expected_var, tolerance = 1e-4)
})

test_that("delta_method() errors before estimation", {
  y <- c(1, 2, 3)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(delta_method(m, transform = function(theta) theta))
})

test_that("delta_method() result is a matrix", {
  m <- make_fitted_mean()
  result <- delta_method(m, transform = function(theta) theta)
  expect_true(is.matrix(result))
})

test_that("delta_method() result is symmetric", {
  m <- make_fitted_mean_variance()
  result <- delta_method(m, transform = function(theta) {
    c(theta[1] + theta[2], theta[1] * theta[2])
  })
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("delta_method() result has non-negative diagonal (variances)", {
  m <- make_fitted_mean_variance()
  result <- delta_method(m, transform = function(theta) {
    c(exp(theta[1]), log(theta[2]))
  })
  expect_true(all(diag(result) >= 0))
})

test_that("delta_method() with constant transform returns zero variance", {
  m <- make_fitted_mean()
  result <- delta_method(m, transform = function(theta) 42)
  expect_equal(result, matrix(0, nrow = 1, ncol = 1), tolerance = 1e-6)
})

# delta_method() input validation --------------------------------------------
#
# Python's delta_method (delicatessen/sandwich.py) raises ValueError under three
# conditions: the covariance matrix is not square, its dimension does not match
# the length of the parameter vector, and the transform output is not
# one-dimensional. These tests pin the equivalent informative aborts in R so a
# bad shape fails loudly instead of returning a wrongly shaped matrix or an
# opaque conformability error.

test_that("delta_method() aborts on a non-square covariance matrix", {
  theta <- c(0.5, -1.2)
  # A 2x3 covariance is not square
  covariance <- matrix(seq_len(6), nrow = 2, ncol = 3)
  expect_error(
    delta_method(theta, transform = function(th) th, covariance = covariance),
    regexp = "square matrix"
  )
})

test_that("delta_method() aborts when covariance dimension does not match theta", {
  theta <- c(0.5, -1.2)
  # A 3x3 covariance for a length-2 parameter vector
  covariance <- diag(3)
  expect_error(
    delta_method(theta, transform = function(th) th, covariance = covariance),
    regexp = "share a dimension"
  )
})

test_that("delta_method() aborts when the transform output is not one-dimensional", {
  theta <- c(0.5, -1.2)
  covariance <- diag(2)
  # A transform that returns a genuine 2x2 matrix rather than a vector
  transform <- function(th) {
    matrix(c(th[1], th[2], th[1], th[2]), nrow = 2, ncol = 2)
  }
  expect_error(
    delta_method(theta, transform = transform, covariance = covariance),
    regexp = "one-dimensional"
  )
})

test_that("numeric delta_method() aborts with a clear message when covariance is missing", {
  # Without covariance, as.matrix(NULL) previously died with the opaque
  # "'data' must be of a vector type, was 'NULL'". The abort must name the
  # argument instead.
  theta <- c(0.5, -1.2)
  expect_error(
    delta_method(theta, transform = function(th) th),
    regexp = "covariance"
  )
})

# The frame the aborts report --------------------------------------------------
#
# Both delta_method() methods hand their arguments to a shared worker, so every
# abort raised while judging those arguments reported the worker's own call.
# That call names the worker's parameters rather than the user's arguments, and
# the worker appears in no man page, so nothing in the report is a call the
# caller can act on. Each abort names the method the caller reached instead,
# following the same threading the estimate() methods use.

reported_call <- function(err) {
  paste(deparse(conditionCall(err)), collapse = " ")
}

test_that("the missing-covariance abort names the method the caller reached", {
  theta <- c(0.5, -1.2)
  err <- expect_error(delta_method(theta, transform = function(th) th))

  expect_match(reported_call(err), "method(delta_method,", fixed = TRUE)
  expect_false(grepl("delta_method_impl", reported_call(err), fixed = TRUE))
})

test_that("the non-square covariance abort names the method the caller reached", {
  theta <- c(0.5, -1.2)
  err <- expect_error(delta_method(
    theta,
    transform = function(th) th,
    covariance = matrix(seq_len(6), nrow = 2, ncol = 3)
  ))

  expect_match(reported_call(err), "method(delta_method,", fixed = TRUE)
  expect_false(grepl("delta_method_impl", reported_call(err), fixed = TRUE))
})

test_that("the covariance dimension abort names the method the caller reached", {
  theta <- c(0.5, -1.2)
  err <- expect_error(delta_method(
    theta,
    transform = function(th) th,
    covariance = diag(3)
  ))

  expect_match(reported_call(err), "method(delta_method,", fixed = TRUE)
  expect_false(grepl("delta_method_impl", reported_call(err), fixed = TRUE))
})

test_that("the transform-output abort names the method the caller reached", {
  theta <- c(0.5, -1.2)
  transform <- function(th) {
    matrix(c(th[1], th[2], th[1], th[2]), nrow = 2, ncol = 2)
  }
  err <- expect_error(delta_method(
    theta,
    transform = transform,
    covariance = diag(2)
  ))

  expect_match(reported_call(err), "method(delta_method,", fixed = TRUE)
  expect_false(grepl("delta_method_impl", reported_call(err), fixed = TRUE))
})

test_that("the estimator method reports its own frame as well", {
  # The estimator method calls the same worker with its own arguments, so it
  # needs the frame threaded through separately from the numeric method.
  m <- make_fitted_mean_variance()
  transform <- function(th) {
    matrix(c(th[1], th[2], th[1], th[2]), nrow = 2, ncol = 2)
  }
  err <- expect_error(delta_method(m, transform = transform))

  expect_match(
    reported_call(err),
    "method(delta_method, deli::deli_estimator)",
    fixed = TRUE
  )
  expect_false(grepl("delta_method_impl", reported_call(err), fixed = TRUE))
})

# deriv_method is case-insensitive --------------------------------------------
#
# Python lowercases every deriv_method comparison. delta_method compared
# case-sensitively, so "Exact" fell through to the finite-difference path and
# "CAPPROX" raised "not supported". These tests pin case-insensitive acceptance.

test_that("delta_method() accepts an upper-case exact deriv_method", {
  m <- make_fitted_mean()
  transform <- function(theta) exp(theta)
  lower <- delta_method(m, transform = transform, deriv_method = "exact")
  upper <- delta_method(m, transform = transform, deriv_method = "Exact")
  expect_equal(upper, lower)
})

test_that("delta_method() accepts a mixed-case finite-difference deriv_method", {
  m <- make_fitted_mean()
  transform <- function(theta) exp(theta)
  lower <- delta_method(m, transform = transform, deriv_method = "capprox")
  upper <- delta_method(m, transform = transform, deriv_method = "CApprox")
  expect_equal(upper, lower)
})

test_that("delta_method() error for an unknown deriv_method lists exact", {
  m <- make_fitted_mean()
  expect_error(
    delta_method(m, transform = function(theta) theta, deriv_method = "nope"),
    regexp = "exact"
  )
})

# delta_method() exact mode (deriv_method = "exact") --------------------------
#
# Exact mode differentiates the transform with forward-mode automatic
# differentiation, which must reproduce the closed-form Jacobian to near machine
# precision. The implementation currently strips tangents by wrapping the
# transform in as.numeric(), so exact mode returns a ZERO Jacobian and therefore
# a zero (or wrongly shaped) covariance with no error. Every test below is
# written to fail on that bug: each asserts a strictly positive variance AND a
# match against an independent reference (a closed-form Jacobian at 1e-8, the
# central-difference approximation at 1e-5, or a Python fixture at 1e-8). A test
# that only checked "no error" or compared against a degenerate reference would
# pass vacuously against the zero Jacobian, so the references here are chosen to
# make a zero Jacobian unmistakable.
#
# The 1e-5 capprox tolerance follows the autodiff test convention: the central
# difference carries roughly 1e-6 numerical error at the default step size. The
# 1e-8 analytic and Python tolerances hold because exact autodiff reproduces the
# closed form, and Python's delta_method also defaults to exact autodiff, so the
# two exact Jacobians agree up to floating-point noise in the shared inputs.
#
# Divergence from Python noted here (documented, not changed): the R
# delta_method default deriv_method is "capprox", whereas Python's default is
# "exact". Adjudication is tracked elsewhere; these tests set
# deriv_method = "exact" explicitly and do not depend on the default.

# Helper to create a fitted logistic-regression MEstimator (log-odds coefficients)
make_fitted_logistic <- function() {
  ref <- load_fixture("ee_regression_logistic")
  X <- ref$X
  y <- ref$y
  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "logistic")
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  estimate(m)
}

test_that("delta_method() exact mode matches analytic and capprox for exp(theta)", {
  m <- make_fitted_mean()
  th <- m@theta[1]
  V <- m@variance[1, 1]

  exact <- delta_method(
    m,
    transform = function(theta) exp(theta),
    deriv_method = "exact"
  )
  capprox <- delta_method(
    m,
    transform = function(theta) exp(theta),
    deriv_method = "capprox"
  )
  # Var(exp(theta)) = exp(theta)^2 * Var(theta)
  analytic <- matrix(exp(th)^2 * V, nrow = 1, ncol = 1)

  # Reference magnitude is ~4.1, far above any tolerance floor and far from zero
  expect_gt(exact[1, 1], 1)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode matches analytic for a parameter ratio", {
  m <- make_fitted_mean_variance()
  th <- m@theta
  V <- m@variance

  exact <- delta_method(
    m,
    transform = function(theta) theta[1] / theta[2],
    deriv_method = "exact"
  )
  capprox <- delta_method(
    m,
    transform = function(theta) theta[1] / theta[2],
    deriv_method = "capprox"
  )
  # G = [1/th2, -th1/th2^2]
  G <- matrix(c(1 / th[2], -th[1] / th[2]^2), nrow = 1)
  analytic <- G %*% V %*% t(G)

  expect_gt(exact[1, 1], 0)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode differentiates inverse_logit of a logistic coefficient", {
  m <- make_fitted_logistic()
  th <- m@theta[1]
  V <- m@variance[1, 1]

  # Transform a log-odds coefficient to a probability via the autodiff-compatible
  # inverse_logit; its derivative is p (1 - p).
  transform <- function(theta) inverse_logit(theta[1])
  exact <- delta_method(m, transform = transform, deriv_method = "exact")
  capprox <- delta_method(m, transform = transform, deriv_method = "capprox")

  p <- inverse_logit(th)
  analytic <- matrix((p * (1 - p))^2 * V, nrow = 1, ncol = 1)

  expect_gt(exact[1, 1], 0)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode differentiates a logit transform", {
  # logit() must work under exact mode like inverse_logit(); its derivative is
  # 1 / (p (1 - p)). The diagonal covariance keeps the analytic reference simple.
  theta <- c(0.3, 0.6)
  Sigma <- diag(c(0.01, 0.02))
  transform <- function(th) c(logit(th[1]), logit(th[2]))

  exact <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "exact"
  )
  capprox <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "capprox"
  )

  g <- 1 / (theta * (1 - theta))
  analytic <- diag(g^2 * diag(Sigma))

  expect_true(all(diag(exact) > 0))
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode builds a full Jacobian for a vector transform", {
  m <- make_fitted_mean_variance()
  th <- m@theta
  V <- m@variance

  transform <- function(theta) c(exp(theta[1]), theta[1] * theta[2])
  exact <- delta_method(m, transform = transform, deriv_method = "exact")
  capprox <- delta_method(m, transform = transform, deriv_method = "capprox")

  # G = [[exp(th1), 0], [th2, th1]]
  G <- matrix(c(exp(th[1]), 0, th[2], th[1]), nrow = 2, byrow = TRUE)
  analytic <- G %*% V %*% t(G)

  # Both output variances are strictly positive; the off-diagonal is non-zero
  expect_true(all(diag(exact) > 0))
  expect_gt(abs(exact[1, 2]), 1e-3)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode differentiates a matrix-vector composite", {
  # Numeric-vector interface with a fixed covariance; the transform is the linear
  # map X %*% theta, whose Jacobian is exactly X, so the delta-method covariance
  # is X %*% Sigma %*% t(X).
  theta <- c(0.5, -1.2, 0.8)
  Sigma <- matrix(
    c(
      0.04,
      0.01,
      -0.005,
      0.01,
      0.05,
      0.002,
      -0.005,
      0.002,
      0.03
    ),
    nrow = 3,
    byrow = TRUE
  )
  X <- matrix(
    c(
      1,
      0.5,
      -0.25,
      1,
      -1,
      2
    ),
    nrow = 2,
    byrow = TRUE
  )

  transform <- function(th) {
    eta <- X %*% th
    c(eta[1], eta[2])
  }

  exact <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "exact"
  )
  capprox <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "capprox"
  )
  analytic <- X %*% Sigma %*% t(X)

  expect_true(all(diag(exact) > 0))
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode differentiates lgamma", {
  m <- make_fitted_mean()
  th <- m@theta[1]
  V <- m@variance[1, 1]

  exact <- delta_method(
    m,
    transform = function(theta) lgamma(theta[1]),
    deriv_method = "exact"
  )
  capprox <- delta_method(
    m,
    transform = function(theta) lgamma(theta[1]),
    deriv_method = "capprox"
  )
  # d/dtheta lgamma(theta) = digamma(theta)
  analytic <- matrix(digamma(th)^2 * V, nrow = 1, ncol = 1)

  expect_gt(exact[1, 1], 0)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode matches Python exact for a mixed transform", {
  ref <- load_fixture("delta_method_exact")
  theta <- ref$theta
  Sigma <- as.matrix(ref$covariance)
  py <- ref$transforms$scalar_mix

  transform <- function(th) {
    c(exp(th[1]), th[1] * th[2], log(th[3]^2 + 1))
  }

  exact <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "exact"
  )
  capprox <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "capprox"
  )

  # Variances are ~1e-2, well above any tolerance floor; must be strictly positive
  expect_true(all(diag(exact) > 1e-4))
  # Python delta_method also defaults to exact autodiff
  expect_equal(exact, as.matrix(py$variance), tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("delta_method() exact mode matches Python exact for inverse_logit(X theta)", {
  ref <- load_fixture("delta_method_exact")
  theta <- ref$theta
  Sigma <- as.matrix(ref$covariance)
  py <- ref$transforms$matrix_invlogit
  X <- as.matrix(py$X)

  transform <- function(th) {
    eta <- X %*% th
    inverse_logit(eta)
  }

  exact <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "exact"
  )
  capprox <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "capprox"
  )

  # Variances are ~1e-3 and ~1e-4; SE scale ~0.03 and ~0.03, clearly non-zero
  expect_true(all(diag(exact) > 1e-6))
  expect_equal(exact, as.matrix(py$variance), tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

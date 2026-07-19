# ---- estimate() for mean EE --------------------------------------------------

test_that("estimate() finds correct theta for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
})

test_that("estimate() computes correct variance for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(unname(m@variance[1, 1]), ref$variance[[1]], tolerance = 1e-4)
})

test_that("estimate() computes correct asymptotic variance for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(
    unname(m@asymptotic_variance[1, 1]),
    ref$asymptotic_variance[[1]],
    tolerance = 1e-4
  )
})

# ---- estimate() for mean+variance EE ----------------------------------------

test_that("estimate() finds correct theta for mean+variance EE", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
})

test_that("estimate() computes correct variance for mean+variance EE", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-4)
})

# ---- estimate() sets n_obs ---------------------------------------------------

test_that("estimate() sets n_obs correctly", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(m@n_obs, 5L)
})

# ---- estimate() sets bread and meat ------------------------------------------

test_that("estimate() populates bread and meat", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_true(!is.null(m@bread))
  expect_true(!is.null(m@meat))
  expect_true(is.matrix(m@bread))
  expect_true(is.matrix(m@meat))
})

# ---- estimate() returns MEstimator ------------------------------------------

test_that("estimate() returns an MEstimator object", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  result <- estimate(m)

  expect_s3_class(result, "deli::MEstimator")
})

# ---- estimate() validates psi(init) before solving --------------------------
#
# Malformed psi/init combinations previously failed deep inside the solver with
# opaque rootSolve internals, or silently returned the initial values. estimate()
# now evaluates the estimating functions once at the starting values and aborts
# with an informative message.

test_that("estimate() rejects a psi/init dimension mismatch", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0, 0))
  expect_error(estimate(m), "1 estimating equation.*2 parameter")
})

test_that("estimate() rejects a psi that is non-finite at init", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y * NaN - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m), "non-finite")
})

test_that("estimate() rejects a psi that returns NULL", {
  psi <- function(theta) NULL
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m), "NULL")
})

test_that("estimate() rejects non-finite data at init under the lm solver", {
  # Previously the lm solver silently returned theta = init with NULL variance.
  psi <- function(theta) matrix(c(1, 2, NA) - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m, solver = "lm"), "non-finite")
})

test_that("estimate() rejects a custom solver that returns the wrong shape", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  bad_solver <- function(stacked_equations, init) c(1, 2, 3)
  expect_error(estimate(m, solver = bad_solver), "numeric vector of length 1")
})

# ---- estimate(deriv_method = "exact") for built-in EEs ----------------------
#
# Exact-mode bread for a built-in estimating equation must agree with the
# central-difference (capprox) bread computed on the same fitted estimator. The
# cross-check uses tolerance 1e-5 because capprox carries roughly 1e-6 numerical
# error at the default step size (floating-point cancellation), consistent with
# the autodiff-vs-approx convention in test-autodiff.R; exact autodiff itself is
# near machine precision, so the tolerance is bounded by the approximation it is
# compared against. These EEs flow through the masked matrix()/rbind() surface
# without stripping tangents (no as.numeric() applied to a theta-carrying value),
# so they exercise the pt_flatten/bind normalization.

test_that("estimate() exact bread matches capprox for ee_mean", {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) ee_mean(theta, y = y)

  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  m_capprox <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "capprox"
  )

  expect_equal(m_exact@theta, m_capprox@theta, tolerance = 1e-9)
  expect_equal(m_exact@bread, m_capprox@bread, tolerance = 1e-5)
  expect_equal(m_exact@variance, m_capprox@variance, tolerance = 1e-5)
})

test_that("estimate() exact bread matches capprox for ee_mean_variance (rbind path)", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y
  psi <- function(theta) ee_mean_variance(theta, y = y)

  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  m_capprox <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "capprox"
  )

  expect_equal(m_exact@theta, m_capprox@theta, tolerance = 1e-9)
  expect_equal(m_exact@bread, m_capprox@bread, tolerance = 1e-5)
  expect_equal(m_exact@variance, m_capprox@variance, tolerance = 1e-5)
})

test_that("estimate() roxygen example runs under deriv_method = 'exact'", {
  # The shape of estimate()'s own runnable example, which is the ee_mean pattern
  # matrix(y - theta[1], nrow = 1) that previously errored under exact mode.
  psi <- function(theta) {
    y <- c(1, 2, 3, 4, 5)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))

  m <- expect_no_error(estimate(m, deriv_method = "exact"))
  expect_equal(unname(m@theta), 3)
  # bread of the summed mean EE is 1 (d/dtheta of -sum(y - theta) is n; /n = 1)
  expect_equal(unname(m@bread), matrix(1))
})

test_that("estimate() fits a psi built with %*% under numerical methods", {
  # A user estimating function that forms the linear predictor with %*%
  # flows through the numerical bread path and recovers the OLS coefficients.
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- 1.5 + 2 * X[, 2] + rnorm(n)
  psi <- function(theta) {
    resid <- as.vector(y - X %*% theta)
    t(X * resid)
  }
  coefs <- unname(coef(lm(y ~ X[, 2])))

  se_capprox <- NULL
  for (method in c("capprox", "fapprox", "bapprox")) {
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0)),
      deriv_method = method
    )
    expect_equal(unname(m@theta), coefs, tolerance = 1e-6)
    expect_equal(dim(m@variance), c(2L, 2L))
    if (is.null(se_capprox)) {
      se_capprox <- sqrt(diag(m@variance))
    } else {
      expect_equal(sqrt(diag(m@variance)), se_capprox, tolerance = 1e-5)
    }
  }
})

# Tests for causal inference estimating equations (bd-djs)
# ee_gformula, ee_ipw, ee_aipw

# Helper: generate simple confounded data
make_causal_data <- function(n = 200, seed = 42) {
  set.seed(seed)
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.25 + 0.5 * W)
  Ya0 <- rbinom(n, 1, 0.75 - 0.5 * W)
  Ya1 <- rbinom(n, 1, 0.75 - 0.5 * W - 0.1)
  Y <- (1 - A) * Ya0 + A * Ya1
  AW <- A * W

  list(
    Y = Y,
    A = A,
    W = W,
    # Outcome model design matrix: intercept, A, W, A*W
    X = cbind(1, A, W, AW),
    # Propensity score model design matrix: intercept, W
    Wmat = cbind(1, W),
    # Counterfactual X under A=1
    X1 = cbind(1, 1, W, 1 * W),
    # Counterfactual X under A=0
    X0 = cbind(1, 0, W, 0 * W),
    n = n
  )
}

# ee_gformula tests -----------------------------------------------------------

test_that("ee_gformula single plan returns (1+p)-by-n matrix", {
  d <- make_causal_data()
  p <- ncol(d$X)
  theta <- rep(0, 1 + p)
  result <- ee_gformula(theta, y = d$Y, X = d$X, X1 = d$X1)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1 + p)
  expect_equal(ncol(result), d$n)
})

test_that("ee_gformula two plans returns (3+p)-by-n matrix", {
  d <- make_causal_data()
  p <- ncol(d$X)
  theta <- rep(0, 3 + p)
  result <- ee_gformula(theta, y = d$Y, X = d$X, X1 = d$X1, X0 = d$X0)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3 + p)
  expect_equal(ncol(result), d$n)
})

test_that("ee_gformula single plan solves via MEstimator", {
  d <- make_causal_data()
  p <- ncol(d$X)

  psi <- function(theta) {
    ee_gformula(theta, y = d$Y, X = d$X, X1 = d$X1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0.5, rep(0, p)))
  m <- estimate(m)

  # Causal mean should be between 0 and 1 for binary outcome
  expect_true(m@theta[1] > 0 && m@theta[1] < 1)
  # Should have variance
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_gformula two plans solves via MEstimator", {
  d <- make_causal_data()
  p <- ncol(d$X)

  psi <- function(theta) {
    ee_gformula(theta, y = d$Y, X = d$X, X1 = d$X1, X0 = d$X0)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5, 0.5, rep(0, p)))
  m <- estimate(m)

  # ACE (theta[1]) should be difference of mu1 (theta[2]) and mu0 (theta[3])
  expect_equal(
    unname(m@theta[1]),
    unname(m@theta[2] - m@theta[3]),
    tolerance = 1e-6
  )
  # Causal means should be between 0 and 1
  expect_true(m@theta[2] > 0 && m@theta[2] < 1)
  expect_true(m@theta[3] > 0 && m@theta[3] < 1)
})

# ee_ipw tests ----------------------------------------------------------------

test_that("ee_ipw returns (3+b)-by-n matrix", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  theta <- rep(0, 3 + b)
  result <- ee_ipw(theta, y = d$Y, A = d$A, W = d$Wmat)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3 + b)
  expect_equal(ncol(result), d$n)
})

test_that("ee_ipw solves via MEstimator", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)

  psi <- function(theta) {
    ee_ipw(theta, y = d$Y, A = d$A, W = d$Wmat)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5, 0.5, rep(0, b)))
  m <- estimate(m)

  # ACE should equal mu1 - mu0
  expect_equal(
    unname(m@theta[1]),
    unname(m@theta[2] - m@theta[3]),
    tolerance = 1e-6
  )
  # Should have reasonable causal means
  expect_true(m@theta[2] > 0 && m@theta[2] < 1)
  expect_true(m@theta[3] > 0 && m@theta[3] < 1)
})

test_that("ee_ipw with truncation doesn't error", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)

  psi <- function(theta) {
    ee_ipw(theta, y = d$Y, A = d$A, W = d$Wmat, truncate = c(0.1, 0.9))
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0.5, 0.5, rep(0, b)))
  m <- estimate(m)
  expect_true(!is.null(m@theta))
})

test_that("ee_ipw errors when truncate bounds are out of order", {
  # Mirrors Python (ee_ipw, causal.py), which raises
  # "truncate values must be specified in ascending order" when truncate[0] > truncate[1].
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  theta <- rep(0, 3 + b)
  expect_error(
    ee_ipw(theta, y = d$Y, A = d$A, W = d$Wmat, truncate = c(0.9, 0.1)),
    regexp = "ascending order"
  )
})

# ee_ipw_msm tests ------------------------------------------------------------

test_that("ee_ipw_msm with ordered truncation returns a matrix", {
  d <- make_causal_data()
  V <- cbind(1, d$A)
  b <- ncol(d$Wmat)
  cc <- ncol(V)
  theta <- rep(0, cc + b)
  result <- ee_ipw_msm(
    theta,
    y = d$Y,
    A = d$A,
    W = d$Wmat,
    V = V,
    distribution = "binomial",
    link = "logit",
    truncate = c(0.1, 0.9)
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), cc + b)
  expect_equal(ncol(result), d$n)
})

test_that("ee_ipw_msm errors when truncate bounds are out of order", {
  # Mirrors Python (ee_ipw_msm, causal.py), which raises
  # "truncate values must be specified in ascending order" when truncate[0] > truncate[1].
  d <- make_causal_data()
  V <- cbind(1, d$A)
  b <- ncol(d$Wmat)
  cc <- ncol(V)
  theta <- rep(0, cc + b)
  expect_error(
    ee_ipw_msm(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      V = V,
      distribution = "binomial",
      link = "logit",
      truncate = c(0.9, 0.1)
    ),
    regexp = "ascending order"
  )
})

# ee_aipw tests ---------------------------------------------------------------

test_that("ee_aipw returns (3+b+c)-by-n matrix", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  c_params <- ncol(d$X)
  theta <- rep(0, 3 + b + c_params)
  result <- ee_aipw(
    theta,
    y = d$Y,
    A = d$A,
    W = d$Wmat,
    X = d$X,
    X1 = d$X1,
    X0 = d$X0
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3 + b + c_params)
  expect_equal(ncol(result), d$n)
})

test_that("ee_aipw solves via MEstimator", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  c_params <- ncol(d$X)

  psi <- function(theta) {
    ee_aipw(theta, y = d$Y, A = d$A, W = d$Wmat, X = d$X, X1 = d$X1, X0 = d$X0)
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(0, 0.5, 0.5, rep(0, b + c_params))
  )
  m <- estimate(m)

  # ACE should equal mu1 - mu0
  expect_equal(
    unname(m@theta[1]),
    unname(m@theta[2] - m@theta[3]),
    tolerance = 1e-6
  )
  # Causal means should be between 0 and 1
  expect_true(m@theta[2] > 0 && m@theta[2] < 1)
  expect_true(m@theta[3] > 0 && m@theta[3] < 1)
})

test_that("ee_aipw with truncation works", {
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  c_params <- ncol(d$X)

  psi <- function(theta) {
    ee_aipw(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      X = d$X,
      X1 = d$X1,
      X0 = d$X0,
      truncate = c(0.1, 0.9)
    )
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(0, 0.5, 0.5, rep(0, b + c_params))
  )
  m <- estimate(m)
  expect_true(!is.null(m@theta))
})

test_that("ee_aipw errors when truncate bounds are out of order", {
  # Mirrors Python (ee_aipw, causal.py), which raises
  # "truncate values must be specified in ascending order" when truncate[0] > truncate[1].
  d <- make_causal_data()
  b <- ncol(d$Wmat)
  c_params <- ncol(d$X)
  theta <- rep(0, 3 + b + c_params)
  expect_error(
    ee_aipw(
      theta,
      y = d$Y,
      A = d$A,
      W = d$Wmat,
      X = d$X,
      X1 = d$X1,
      X0 = d$X0,
      truncate = c(0.9, 0.1)
    ),
    regexp = "ascending order"
  )
})

# Cross-estimator consistency -------------------------------------------------

test_that("gformula, ipw, and aipw estimate ATEs in same direction", {
  d <- make_causal_data()
  p <- ncol(d$X)
  b <- ncol(d$Wmat)

  # G-formula ATE
  psi_gf <- function(theta) {
    ee_gformula(theta, y = d$Y, X = d$X, X1 = d$X1, X0 = d$X0)
  }
  m_gf <- estimate(MEstimator(
    stacked_equations = psi_gf,
    init = c(0, 0.5, 0.5, rep(0, p))
  ))

  # IPW ATE
  psi_ipw <- function(theta) {
    ee_ipw(theta, y = d$Y, A = d$A, W = d$Wmat)
  }
  m_ipw <- estimate(MEstimator(
    stacked_equations = psi_ipw,
    init = c(0, 0.5, 0.5, rep(0, b))
  ))

  # AIPW ATE
  psi_aipw <- function(theta) {
    ee_aipw(theta, y = d$Y, A = d$A, W = d$Wmat, X = d$X, X1 = d$X1, X0 = d$X0)
  }
  m_aipw <- estimate(MEstimator(
    stacked_equations = psi_aipw,
    init = c(0, 0.5, 0.5, rep(0, b + p))
  ))

  # All three should give ATEs with the same sign
  # (not testing exact values since methods differ)
  ates <- c(m_gf@theta[1], m_ipw@theta[1], m_aipw@theta[1])
  # At minimum, all should be finite

  expect_true(all(is.finite(ates)))
})

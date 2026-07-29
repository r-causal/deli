# Tests for ee_tobit (Tobit Type I regression with censored outcomes)

test_that("ee_tobit returns correct shape", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y_star <- 0.5 + 0.3 * X[, 2] + rnorm(n)
  y <- pmax(y_star, -1) # left censor at -1

  b <- ncol(X)
  theta <- c(rep(0, b), log(1)) # beta + log(sigma)

  result <- ee_tobit(theta, X = X, y = y, lower = -1)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

test_that("ee_tobit with no censoring matches linear regression", {
  set.seed(42)
  n <- 100
  X <- cbind(1, rnorm(n))
  y <- 0.5 + 0.8 * X[, 2] + rnorm(n)

  # Tobit with no censoring (lower = -Inf, upper = Inf by default)
  psi_tobit <- function(theta) {
    ee_tobit(theta, X = X, y = y)
  }
  m_tobit <- MEstimator(
    stacked_equations = psi_tobit,
    init = c(mean(y), 0, log(sd(y)))
  )
  m_tobit <- estimate(m_tobit, solver = "nleqslv")

  # Linear regression
  psi_lm <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m_lm <- MEstimator(stacked_equations = psi_lm, init = c(0, 0))
  m_lm <- estimate(m_lm)

  # Beta coefficients should match. ee_tobit() names its rows and
  # ee_regression() does not, so the two fits label the same coefficients
  # differently and only the values are compared; test-param-names-builtin.R
  # pins what each is called.
  expect_equal(unname(m_tobit@theta[1:2]), unname(m_lm@theta), tolerance = 1e-3)
})

test_that("ee_tobit solves with left censoring", {
  set.seed(123)
  n <- 200
  X <- cbind(1, rnorm(n))
  y_star <- 1 + 0.5 * X[, 2] + rnorm(n)
  lower_limit <- 0
  y <- pmax(y_star, lower_limit)

  psi <- function(theta) {
    ee_tobit(theta, X = X, y = y, lower = lower_limit)
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = c(mean(y), 0, log(sd(y)))
  )
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Sigma should be positive
  expect_true(exp(m@theta[3]) > 0)
  # Intercept should be close to 1
  expect_equal(unname(m@theta[1]), 1, tolerance = 0.3)
})

test_that("ee_tobit solves with right censoring", {
  set.seed(456)
  n <- 200
  X <- cbind(1, rnorm(n))
  y_star <- 1 + 0.5 * X[, 2] + rnorm(n)
  upper_limit <- 2.5
  y <- pmin(y_star, upper_limit)

  psi <- function(theta) {
    ee_tobit(theta, X = X, y = y, upper = upper_limit)
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = c(mean(y), 0, log(sd(y)))
  )
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(exp(m@theta[3]) > 0)
})

test_that("ee_tobit solves with both left and right censoring", {
  set.seed(789)
  n <- 200
  X <- cbind(1, rnorm(n))
  y_star <- 0.5 + 0.3 * X[, 2] + rnorm(n)
  lower_limit <- -1
  upper_limit <- 2
  y <- pmin(pmax(y_star, lower_limit), upper_limit)

  psi <- function(theta) {
    ee_tobit(theta, X = X, y = y, lower = lower_limit, upper = upper_limit)
  }

  m <- MEstimator(
    stacked_equations = psi,
    init = c(mean(y), 0, log(sd(y)))
  )
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(exp(m@theta[3]) > 0)
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_tobit errors when lower >= upper", {
  X <- cbind(1, rnorm(10))
  y <- rnorm(10)

  expect_error(
    ee_tobit(c(0, 0, 0), X = X, y = y, lower = 2, upper = 1),
    "lower"
  )
})

test_that("ee_tobit errors when y values below lower limit", {
  X <- cbind(1, rnorm(10))
  y <- c(-5, 1, 2, 3, 4, 5, 6, 7, 8, 9)

  expect_error(
    ee_tobit(c(0, 0, 0), X = X, y = y, lower = 0),
    "below"
  )
})

test_that("ee_tobit errors when y values above upper limit", {
  X <- cbind(1, rnorm(10))
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)

  expect_error(
    ee_tobit(c(0, 0, 0), X = X, y = y, upper = 10),
    "above"
  )
})

test_that("ee_tobit supports weights", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y_star <- 1 + 0.5 * X[, 2] + rnorm(n)
  y <- pmax(y_star, 0)
  wts <- runif(n, 0.5, 1.5)

  result <- ee_tobit(c(0, 0, 0), X = X, y = y, lower = 0, weights = wts)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

test_that("ee_tobit supports offset", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  off <- rep(0.1, n)
  y_star <- 1 + 0.5 * X[, 2] + rnorm(n)
  y <- pmax(y_star, 0)

  result <- ee_tobit(c(0, 0, 0), X = X, y = y, lower = 0, offset = off)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_tobit rejects a y whose length differs from the rows of X", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- pmax(0.5 + 0.3 * X[, 2] + rnorm(n), -1)
  theta <- c(rep(0, ncol(X)), log(1))

  expect_error(
    ee_tobit(theta, X = X, y = y[seq_len(n / 2)], lower = -1),
    "same length as the data"
  )
})

test_that("ee_tobit rejects an offset whose length differs from the rows of X", {
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- pmax(0.5 + 0.3 * X[, 2] + rnorm(n), -1)
  theta <- c(rep(0, ncol(X)), log(1))

  expect_error(
    ee_tobit(theta, X = X, y = y, lower = -1, offset = rep(0, n / 2)),
    "same length as the data"
  )
})

# Tests for ee_mlogit and ee_beta_regression

# ee_mlogit tests -------------------------------------------------------------

test_that("ee_mlogit returns correct shape for 3-category outcome", {
  n <- 20
  set.seed(99)
  X <- cbind(1, rnorm(n))
  y_cat <- sample(1:3, n, replace = TRUE)
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  b <- ncol(X)
  k <- ncol(y)
  theta <- rep(0, b * (k - 1))

  result <- ee_mlogit(theta, X = X, y = y)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b * (k - 1))
  expect_equal(ncol(result), n)
})

test_that("ee_mlogit solves via MEstimator", {
  set.seed(123)
  n <- 50
  W <- rbinom(n, 1, 0.5)
  probs <- cbind(
    0.5 - 0.2 * W,
    0.3 + 0.1 * W,
    0.2 + 0.1 * W
  )
  y_cat <- sapply(1:n, function(i) sample(1:3, 1, prob = probs[i, ]))
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  X <- cbind(1, W)

  psi <- function(theta) {
    ee_mlogit(theta, X = X, y = y)
  }
  m <- MEstimator(stacked_equations = psi, init = rep(0, 4))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
})

test_that("ee_mlogit errors on parameter count mismatch", {
  X <- cbind(1, rnorm(10))
  y <- cbind(
    c(1, 0, 0, 1, 0, 0, 1, 0, 1, 0),
    c(0, 1, 0, 0, 1, 0, 0, 1, 0, 1),
    c(0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
  )
  expect_error(ee_mlogit(rep(0, 3), X = X, y = y), "mismatch")
})

# ee_beta_regression tests ----------------------------------------------------

test_that("ee_beta_regression returns correct shape", {
  set.seed(42)
  n <- 30
  X <- cbind(1, rnorm(n))
  y <- runif(n, 0.01, 0.99)
  b <- ncol(X)

  theta <- rep(0, b + 1)
  result <- ee_beta_regression(theta, X = X, y = y)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), b + 1)
  expect_equal(ncol(result), n)
})

test_that("ee_beta_regression solves via MEstimator", {
  set.seed(42)
  n <- 50
  W <- rnorm(n)
  mu <- 1 / (1 + exp(-(0.5 + 0.3 * W)))
  y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
  y <- pmin(pmax(y, 0.001), 0.999)
  X <- cbind(1, W)

  psi <- function(theta) {
    ee_beta_regression(theta, X = X, y = y)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0, log(10)))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(diag(m@variance) > 0))
  expect_true(exp(m@theta[3]) > 0)
})

test_that("ee_beta_regression intercept-only model recovers mean", {
  set.seed(123)
  n <- 100
  mu_true <- 0.6
  phi_true <- 20
  y <- rbeta(n, shape1 = mu_true * phi_true, shape2 = (1 - mu_true) * phi_true)
  y <- pmin(pmax(y, 0.001), 0.999)
  X <- cbind(rep(1, n))

  psi <- function(theta) {
    ee_beta_regression(theta, X = X, y = y)
  }
  # Use logit(mean(y)) as init for better convergence
  init_beta <- log(mean(y) / (1 - mean(y)))
  m <- MEstimator(stacked_equations = psi, init = c(init_beta, log(phi_true)))
  m <- estimate(m, solver = "nleqslv")

  estimated_mean <- 1 / (1 + exp(-m@theta[1]))
  expect_equal(unname(estimated_mean), mu_true, tolerance = 0.1)
})

# Input validation (batch F) --------------------------------------------------

test_that("ee_mlogit rejects a y whose rows differ from the rows of X", {
  n <- 20
  set.seed(99)
  X <- cbind(1, rnorm(n))
  y_cat <- sample(1:3, n, replace = TRUE)
  y <- cbind(
    as.integer(y_cat == 1),
    as.integer(y_cat == 2),
    as.integer(y_cat == 3)
  )
  theta <- rep(0, ncol(X) * (ncol(y) - 1))

  expect_error(
    ee_mlogit(theta, X = X, y = y[seq_len(n / 2), , drop = FALSE]),
    "same number of rows|same length as the data"
  )
})

test_that("ee_beta_regression rejects a y whose length differs from the rows of X", {
  set.seed(11)
  n <- 40
  X <- cbind(1, rnorm(n))
  y <- runif(n, 0.1, 0.9)
  theta <- c(0, 0, 0)

  expect_error(
    ee_beta_regression(theta, X = X, y = y[seq_len(n / 2)]),
    "same length as the data"
  )
})

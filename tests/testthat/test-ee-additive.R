# Tests for ee_additive_regression()

test_that("ee_additive_regression returns correct dimensions", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.5)
  X <- cbind(1, x)
  knots <- c(-1, 0, 1)
  specs <- list(NULL, list(knots = knots, penalty = 5))

  # Expected number of params: 2 original + 2 restricted spline terms = 4
  theta <- rep(0, 4)
  result <- ee_additive_regression(theta, X, y, specs, model = "linear")

  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), n)
})

test_that("ee_additive_regression with no splines matches ee_regression", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n, sd = 0.5)
  X <- cbind(1, x)
  specs <- list(NULL, NULL)

  theta <- c(1, 2)

  result_additive <- ee_additive_regression(
    theta, X, y, specs, model = "linear"
  )
  result_regression <- ee_regression(theta, X, y, model = "linear")

  expect_equal(unname(result_additive), unname(result_regression))
})

test_that("ee_additive_regression with zero penalty on splines matches bridge", {

  set.seed(99)
  n <- 50
  x <- rnorm(n)
  y <- sin(x) + rnorm(n, sd = 0.3)
  X <- cbind(1, x)
  knots <- c(-1, 0, 1)
  specs <- list(NULL, list(knots = knots, penalty = 0))

  # Build additive design matrix manually
  dm <- additive_design_matrix(X, specs, return_penalty = TRUE)
  theta <- rep(0, ncol(dm$X))

  result_additive <- ee_additive_regression(
    theta, X, y, specs, model = "linear"
  )

  result_bridge <- ee_bridge_regression(
    theta, X = dm$X, y = y, model = "linear",
    penalty = dm$penalty, gamma = 2, center = 0
  )

  expect_equal(result_additive, result_bridge)
})

test_that("ee_additive_regression linear GAM can be estimated", {
  set.seed(42)
  n <- 200
  x <- runif(n, -3, 3)
  # Non-linear relationship
  y <- sin(x) + rnorm(n, sd = 0.3)
  X <- cbind(1, x)
  knots <- seq(-2, 2, by = 1)
  specs <- list(NULL, list(knots = knots, penalty = 10))

  dm <- additive_design_matrix(X, specs)
  n_params <- ncol(dm)
  init <- rep(0, n_params)

  psi <- function(theta) {
    ee_additive_regression(theta, X, y, specs, model = "linear")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  # The estimator should converge
  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
  expect_true(all(diag(m@variance) >= 0))
})

test_that("ee_additive_regression logistic GAM can be estimated", {
  set.seed(123)
  n <- 300
  x <- runif(n, -3, 3)
  prob <- 1 / (1 + exp(-(sin(x) + 0.5)))
  y <- rbinom(n, 1, prob)
  X <- cbind(1, x)
  knots <- seq(-2, 2, by = 1)
  specs <- list(NULL, list(knots = knots, penalty = 10))

  dm <- additive_design_matrix(X, specs)
  n_params <- ncol(dm)
  init <- rep(0, n_params)

  psi <- function(theta) {
    ee_additive_regression(theta, X, y, specs, model = "logistic")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
})

test_that("ee_additive_regression poisson GAM can be estimated", {
  set.seed(321)
  n <- 200
  x <- runif(n, -2, 2)
  lambda <- exp(0.5 + 0.5 * sin(x))
  y <- rpois(n, lambda)
  X <- cbind(1, x)
  knots <- c(-1, 0, 1)
  specs <- list(NULL, list(knots = knots, penalty = 5))

  dm <- additive_design_matrix(X, specs)
  n_params <- ncol(dm)
  init <- rep(0, n_params)

  psi <- function(theta) {
    ee_additive_regression(theta, X, y, specs, model = "poisson")
  }

  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(diag(m@variance))))
})

test_that("ee_additive_regression with weights works", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.5)
  X <- cbind(1, x)
  w <- runif(n, 0.5, 1.5)
  knots <- c(-1, 0, 1)
  specs <- list(NULL, list(knots = knots, penalty = 5))

  dm <- additive_design_matrix(X, specs)
  theta <- rep(0, ncol(dm))

  result <- ee_additive_regression(
    theta, X, y, specs, model = "linear", weights = w
  )

  expect_equal(nrow(result), length(theta))
  expect_equal(ncol(result), n)
})

test_that("ee_additive_regression with offset works", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  offset <- rep(0.5, n)
  y <- x^2 + 0.5 + rnorm(n, sd = 0.5)
  X <- cbind(1, x)
  knots <- c(-1, 0, 1)
  specs <- list(NULL, list(knots = knots, penalty = 5))

  dm <- additive_design_matrix(X, specs)
  theta <- rep(0, ncol(dm))

  result <- ee_additive_regression(
    theta, X, y, specs, model = "linear", offset = offset
  )

  expect_equal(nrow(result), length(theta))
  expect_equal(ncol(result), n)
})

test_that("ee_additive_regression penalty only affects spline terms", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- x + rnorm(n, sd = 0.1)
  X <- cbind(1, x)
  knots <- c(-1, 0, 1)

  # High penalty vs zero penalty
  specs_pen <- list(NULL, list(knots = knots, penalty = 100))
  specs_no <- list(NULL, list(knots = knots, penalty = 0))

  theta <- rep(0.5, 4)

  result_pen <- ee_additive_regression(
    theta, X, y, specs_pen, model = "linear"
  )
  result_no <- ee_additive_regression(
    theta, X, y, specs_no, model = "linear"
  )

  # First two rows (intercept and linear term) should be same
  # because their penalty is always 0
  expect_equal(result_pen[1, ], result_no[1, ])
  expect_equal(result_pen[2, ], result_no[2, ])

  # Last two rows (spline terms) should differ because of penalty
  expect_false(isTRUE(all.equal(result_pen[3, ], result_no[3, ])))
})

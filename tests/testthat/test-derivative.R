# ---- linear function ---------------------------------------------------------

test_that("approx_differentiation computes Jacobian for linear function", {
  # f(x) = c(2*x[1] + 3*x[2], x[1] - x[2])
  # Jacobian: [[2, 3], [1, -1]] (constant, exact for all methods)
  f <- function(x) c(2 * x[1] + 3 * x[2], x[1] - x[2])
  theta <- c(1, 1)
  expected <- matrix(c(2, 1, 3, -1), nrow = 2, ncol = 2)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(result, expected, tolerance = 1e-5)
  }
})

# ---- quadratic function -----------------------------------------------------

test_that("approx_differentiation computes Jacobian for quadratic function", {
  # f(x) = c(x[1]^2, x[2]^2)
  # At theta = c(3, 4): Jacobian = [[6, 0], [0, 8]]
  f <- function(x) c(x[1]^2, x[2]^2)
  theta <- c(3, 4)
  expected <- matrix(c(6, 0, 0, 8), nrow = 2, ncol = 2)

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result, expected, tolerance = 1e-5)
})

test_that("quadratic Jacobian has zero off-diagonal elements", {
  f <- function(x) c(x[1]^2, x[2]^2)
  theta <- c(3, 4)

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 2], 0, tolerance = 1e-5)
  expect_equal(result[2, 1], 0, tolerance = 1e-5)
})

# ---- scalar function ---------------------------------------------------------

test_that("approx_differentiation handles scalar function f(x) = x^3", {
  # f(x) = x^3, derivative = 3x^2 = 12 at x = 2
  f <- function(x) x^3
  theta <- 2

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 1], 12, tolerance = 1e-5)
})

# ---- method agreement --------------------------------------------------------

test_that("central, forward, and backward methods agree for smooth functions", {
  f <- function(x) c(sin(x[1]) + x[2]^2, x[1] * x[2])
  theta <- c(1.5, 2.5)

  result_central  <- approx_differentiation(func = f, theta = theta, method = "capprox")
  result_forward  <- approx_differentiation(func = f, theta = theta, method = "fapprox")
  result_backward <- approx_differentiation(func = f, theta = theta, method = "bapprox")

  expect_equal(result_central, result_forward, tolerance = 1e-5)
  expect_equal(result_central, result_backward, tolerance = 1e-5)
})

# ---- dx parameter ------------------------------------------------------------

test_that("smaller dx gives more accurate results for smooth functions", {
  # f(x) = x^2, derivative at x = 5 is 10
  f <- function(x) x^2
  theta <- 5
  exact <- 10

  result_large <- approx_differentiation(func = f, theta = theta, method = "fapprox", dx = 1e-3)
  result_small <- approx_differentiation(func = f, theta = theta, method = "fapprox", dx = 1e-9)

  error_large <- abs(result_large[1, 1] - exact)
  error_small <- abs(result_small[1, 1] - exact)
  expect_true(error_small < error_large)
})

# ---- exponential function ----------------------------------------------------

test_that("approx_differentiation computes derivative of exp(x) at x = 1", {
  f <- function(x) exp(x)
  theta <- 1

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 1], exp(1), tolerance = 1e-5)
})

# ---- column-matrix-returning function ----------------------------------------

test_that("approx_differentiation handles a func that returns a column matrix", {
  # f(th) = X %*% th returns an n-by-1 column matrix in R (any %*% does).
  # The Jacobian of a linear map is X itself, an n-by-p matrix.
  X <- matrix(c(0.3, -1.1, 0.7, 2.0, -0.5, 1.4), nrow = 3, ncol = 2)
  f <- function(th) X %*% th
  theta <- c(0.5, -0.3)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(3L, 2L))
    expect_equal(unname(result), X, tolerance = 1e-5)
  }
})

test_that("approx_differentiation returns n-by-p for vector-valued func", {
  # A plain vector return of length n with p parameters yields an n-by-p matrix.
  f <- function(x) c(2 * x[1] + 3 * x[2], x[1] - x[2], x[1])
  theta <- c(1, 1)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(3L, 2L))
  }
})

test_that("approx_differentiation returns 1-by-1 for scalar func of one parameter", {
  f <- function(x) x^3
  theta <- 2

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(1L, 1L))
    expect_equal(result[1, 1], 12, tolerance = 1e-5)
  }
})

test_that("approx_differentiation handles a length-1 vector return", {
  # A length-1 vector output with two parameters yields a 1-by-2 gradient row.
  f <- function(x) x[1] * x[2]
  theta <- c(3, 5)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(1L, 2L))
    expect_equal(result[1, 1], 5, tolerance = 1e-5)
    expect_equal(result[1, 2], 3, tolerance = 1e-5)
  }
})

# ---- all methods run without error -------------------------------------------

test_that("all three methods run without error", {
  f <- function(x) c(x[1] + x[2], x[1] * x[2])
  theta <- c(2, 3)

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_no_error(
      approx_differentiation(func = f, theta = theta, method = method)
    )
  }
})

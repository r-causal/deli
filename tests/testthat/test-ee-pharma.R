# Tests for pharmacokinetics estimating equations

# ee_emax tests ---------------------------------------------------------------

test_that("ee_emax returns 3-by-n matrix", {
  dose <- c(0, 1, 2, 5, 10, 20, 50, 100)
  response <- c(0, 5, 8, 15, 20, 23, 24.5, 25)
  theta <- c(0, 25, 10)

  result <- ee_emax(theta, dose = dose, response = response)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), length(dose))
})

test_that("ee_emax solves via MEstimator", {
  set.seed(42)
  n <- 50
  dose <- runif(n, 0, 100)
  e0_true <- 1
  emax_true <- 25
  e50_true <- 10
  response <- e0_true +
    (emax_true * dose) / (e50_true + dose) +
    rnorm(n, sd = 1)

  psi <- function(theta) {
    ee_emax(theta, dose = dose, response = response)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 20, 5))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # e0 should be close to 1
  expect_equal(unname(m@theta[1]), e0_true, tolerance = 1)
  # emax should be close to 25
  expect_equal(unname(m@theta[2]), emax_true, tolerance = 5)
  # e50 should be positive
  expect_true(m@theta[3] > 0)
})

test_that("ee_emax at true parameters has near-zero sum", {
  dose <- c(0, 1, 5, 10, 50, 100)
  e0 <- 2
  emax <- 30
  e50 <- 15
  response <- e0 + (emax * dose) / (e50 + dose)

  # At true parameters, residuals should be zero
  result <- ee_emax(c(e0, emax, e50), dose = dose, response = response)
  expect_equal(
    rowSums(result),
    c(e0 = 0, emax = 0, ed50 = 0),
    tolerance = 1e-10
  )
})

test_that("ee_emax with robust loss works", {
  set.seed(42)
  dose <- c(0, 1, 5, 10, 50, 100)
  response <- c(2, 5, 12, 18, 28, 30)

  result <- ee_emax(
    c(2, 30, 15),
    dose = dose,
    response = response,
    loss = "huber",
    k = 1.345
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
})

# ee_emax_ed tests ------------------------------------------------------------

test_that("ee_emax_ed returns 1-by-n matrix", {
  dose <- c(0, 1, 5, 10, 50, 100)
  result <- ee_emax_ed(theta = 15, dose = dose, delta = 0.5, ed50 = 15)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(dose))
})

test_that("ee_emax_ed at delta=0.5 equals ed50", {
  dose <- 1:10
  # ED50 by definition: delta=0.5 -> theta = 0.5/(1-0.5) * ed50 = ed50
  result <- ee_emax_ed(theta = 10, dose = dose, delta = 0.5, ed50 = 10)
  expect_equal(as.numeric(result[1, 1]), 0, tolerance = 1e-10)
})

# ee_loglogistic tests --------------------------------------------------------

test_that("ee_loglogistic returns 4-by-n matrix", {
  dose <- c(0.01, 0.1, 1, 5, 10, 50, 100)
  response <- c(5, 5.5, 8, 15, 20, 24, 25)
  theta <- c(5, 25, 10, 1)

  result <- ee_loglogistic(theta, dose = dose, response = response)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), length(dose))
})

test_that("ee_loglogistic handles zero dose", {
  dose <- c(0, 1, 5, 10, 50)
  response <- c(5, 8, 15, 20, 24)
  theta <- c(5, 25, 10, 1)

  result <- ee_loglogistic(theta, dose = dose, response = response)
  expect_true(all(is.finite(result)))
})

test_that("ee_loglogistic solves via MEstimator", {
  set.seed(42)
  n <- 40
  dose <- exp(runif(n, log(0.1), log(100)))
  e0 <- 5
  emax <- 25
  e50 <- 10
  steep <- 2
  rho <- (dose / e50)^steep
  response <- e0 + (emax - e0) / (1 + rho) + rnorm(n, sd = 0.5)

  psi <- function(theta) {
    ee_loglogistic(theta, dose = dose, response = response)
  }

  m <- MEstimator(stacked_equations = psi, init = c(4, 26, 8, 1.5))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Lower limit should be close to 5
  expect_equal(unname(m@theta[1]), e0, tolerance = 2)
  # Upper limit should be close to 25
  expect_equal(unname(m@theta[2]), emax, tolerance = 3)
})

test_that("ee_loglogistic with robust loss works", {
  dose <- c(0.1, 1, 5, 10, 50, 100)
  response <- c(5, 8, 15, 20, 24, 25)

  result <- ee_loglogistic(
    c(5, 25, 10, 1),
    dose = dose,
    response = response,
    loss = "huber",
    k = 1.345
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
})

# ee_loglogistic_ed tests -----------------------------------------------------

test_that("ee_loglogistic_ed returns 1-by-n matrix", {
  dose <- c(0.1, 1, 5, 10, 50)
  result <- ee_loglogistic_ed(
    theta = 10,
    dose = dose,
    delta = 0.5,
    lower = 5,
    upper = 25,
    ed50 = 10,
    steepness = 1
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(dose))
})

test_that("ee_loglogistic_ed at ed50 with delta=0.5 is near zero", {
  dose <- 1:10
  # At ED50, the response is midway between lower and upper
  result <- ee_loglogistic_ed(
    theta = 10,
    dose = dose,
    delta = 0.5,
    lower = 0,
    upper = 100,
    ed50 = 10,
    steepness = 1
  )
  expect_equal(as.numeric(result[1, 1]), 0, tolerance = 1e-6)
})

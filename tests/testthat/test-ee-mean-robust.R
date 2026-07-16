# Tests for ee_mean_robust() (bd-ehv.3)

test_that("ee_mean_robust with Huber loss matches Python fixture", {
  ref <- load_fixture("ee_mean_robust_huber")
  y <- ref$y
  k <- ref$k

  psi <- function(theta) {
    ee_mean_robust(theta, y = y, k = k, loss = "huber")
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  expect_equal(unname(m@variance), ref$variance, tolerance = 1e-4)
})

test_that("ee_mean_robust with default loss produces same result", {
  ref <- load_fixture("ee_mean_robust_huber")
  y <- ref$y
  k <- ref$k

  # Default loss is "huber", so omitting the argument should match
  psi <- function(theta) {
    ee_mean_robust(theta, y = y, k = k)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  expect_equal(unname(m@variance), ref$variance, tolerance = 1e-4)
})

test_that("ee_mean_robust returns 1-by-n matrix", {
  y <- c(1, 2, 3, 4, 5)
  result <- ee_mean_robust(theta = c(3), y = y, k = 3, loss = "huber")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(y))
})

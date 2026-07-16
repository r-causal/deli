# Tests for influence_functions() S7 generic

make_fitted_mean <- function() {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  estimate(m)
}

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

test_that("influence_functions() returns n-by-p matrix", {
  m <- make_fitted_mean()
  IF <- influence_functions(m)
  expect_true(is.matrix(IF))
  expect_equal(nrow(IF), m@n_obs)
  expect_equal(ncol(IF), length(m@theta))
})

test_that("influence_functions() returns n-by-p matrix for multi-param", {
  m <- make_fitted_mean_variance()
  IF <- influence_functions(m)
  expect_true(is.matrix(IF))
  expect_equal(nrow(IF), m@n_obs)
  expect_equal(ncol(IF), length(m@theta))
})

test_that("influence_functions() columns sum to approximately zero", {
  m <- make_fitted_mean()
  IF <- influence_functions(m)
  # IF should sum to zero when EE sum to zero at theta-hat
  expect_equal(sum(IF), 0, tolerance = 1e-5)
})

test_that("influence_functions() errors before estimation", {
  y <- c(1, 2, 3)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(influence_functions(m))
})

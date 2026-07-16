# ---- Subset solving only solves specified parameters -------------------------

test_that("subset solving only solves for specified parameters", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, only solve for theta[2] (variance)
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # theta[2] should be solved to the reference variance

  expect_equal(unname(m@theta[2]), ref$theta[2], tolerance = 1e-4)
})

# ---- Non-subset parameters remain at init values ----------------------------

test_that("non-subset parameters remain at init values", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, only solve for theta[2]
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # theta[1] should remain at the init value (the true mean)
  expect_equal(unname(m@theta[1]), true_mean)
})

# ---- Results match when subset includes all parameters -----------------------

test_that("results match when subset includes all parameters", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Estimate without subset
  m_full <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m_full <- estimate(m_full)

  # Estimate with subset = all parameters
  m_subset <- MEstimator(
    stacked_equations = psi,
    init = c(0, 1),
    subset = c(1L, 2L)
  )
  m_subset <- estimate(m_subset)

  expect_equal(m_subset@theta, m_full@theta, tolerance = 1e-6)
})

# ---- Subset with mean+variance EE -------------------------------------------

test_that("subset with mean+variance EE matches reference values", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, solve only for variance (theta[2])
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # The full theta should combine the fixed mean and the solved variance
  expect_equal(unname(m@theta[1]), true_mean)
  expect_equal(unname(m@theta[2]), ref$theta[2], tolerance = 1e-4)

  # bread, meat, and variance should still be computed
  expect_true(!is.null(m@bread))
  expect_true(!is.null(m@meat))
  expect_true(!is.null(m@variance))
})

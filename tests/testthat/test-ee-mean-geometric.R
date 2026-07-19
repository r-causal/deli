# Tests for ee_mean_geometric()

test_that("ee_mean_geometric with log_theta=FALSE solves to mean(log(y))", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  # With log_theta=FALSE, theta is the log of the geometric mean
  psi <- function(theta) {
    ee_mean_geometric(theta, y = y, log_theta = FALSE)
  }

  m <- MEstimator(stacked_equations = psi, init = c(1))
  m <- estimate(m)

  expected_log_gm <- mean(log(y))
  expect_equal(unname(m@theta[1]), expected_log_gm, tolerance = 1e-4)
})

test_that("ee_mean_geometric with log_theta=TRUE solves to geometric mean", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  # With log_theta=TRUE (default), theta is the geometric mean itself
  psi <- function(theta) {
    ee_mean_geometric(theta, y = y, log_theta = TRUE)
  }

  expected_gm <- exp(mean(log(y)))
  m <- MEstimator(stacked_equations = psi, init = c(expected_gm))
  m <- estimate(m)

  expect_equal(unname(m@theta[1]), expected_gm, tolerance = 1e-4)
})

test_that("ee_mean_geometric with weights works correctly", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  w <- rep(1, length(y))

  # Uniform weights should match unweighted result
  psi_weighted <- function(theta) {
    ee_mean_geometric(theta, y = y, weights = w, log_theta = FALSE)
  }

  psi_unweighted <- function(theta) {
    ee_mean_geometric(theta, y = y, log_theta = FALSE)
  }

  m_w <- MEstimator(stacked_equations = psi_weighted, init = c(1))
  m_w <- estimate(m_w)

  m_u <- MEstimator(stacked_equations = psi_unweighted, init = c(1))
  m_u <- estimate(m_u)

  expect_equal(m_w@theta, m_u@theta, tolerance = 1e-4)
})

test_that("ee_mean_geometric returns 1-by-n matrix", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  result <- ee_mean_geometric(theta = c(3), y = y, log_theta = FALSE)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(y))
})

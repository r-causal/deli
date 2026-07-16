# check_alpha -------------------------------------------------------------

test_that("check_alpha accepts valid alpha values", {
  expect_no_error(check_alpha(0.05))
  expect_no_error(check_alpha(0.5))
  expect_no_error(check_alpha(0.99))
  expect_no_error(check_alpha(0.01))
  expect_no_error(check_alpha(0.001))
})

test_that("check_alpha rejects alpha outside (0, 1)", {
  expect_error(check_alpha(0), "alpha.*0.*1")
  expect_error(check_alpha(1), "alpha.*0.*1")
  expect_error(check_alpha(-0.5), "alpha.*0.*1")
  expect_error(check_alpha(1.5), "alpha.*0.*1")
  expect_error(check_alpha(-1), "alpha.*0.*1")
  expect_error(check_alpha(2), "alpha.*0.*1")
})

test_that("check_alpha rejects boundary values exactly at 0 and 1", {
  expect_error(check_alpha(0), "alpha")
  expect_error(check_alpha(1), "alpha")
})

# check_estimated ---------------------------------------------------------

test_that("check_estimated accepts a fitted estimator", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  m <- estimate(MEstimator(stacked_equations = psi, init = c(0)))
  expect_no_error(check_estimated(m))
})

test_that("check_estimated rejects an estimator that has not been fitted", {
  psi <- function(theta) matrix(c(1, 2, 3, 4, 5) - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(check_estimated(m), "before calling")
})

# check_penalty_shape -----------------------------------------------------

test_that("check_penalty_shape accepts scalar penalty and center", {

  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(theta, penalty = 0.5, center = 0))
})

test_that("check_penalty_shape accepts penalty and center matching theta length", {
  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(theta, penalty = c(0.1, 0.2, 0.3), center = c(0, 0, 0)))
})

test_that("check_penalty_shape accepts mixed scalar and vector arguments", {
  theta <- c(1, 2, 3)
  # scalar penalty, vector center

  expect_no_error(check_penalty_shape(theta, penalty = 0.5, center = c(0, 0, 0)))
  # vector penalty, scalar center
  expect_no_error(check_penalty_shape(theta, penalty = c(0.1, 0.2, 0.3), center = 0))
})

test_that("check_penalty_shape rejects penalty with wrong length", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, 0.2), center = 0),
    "penalty.*single number.*same length"
  )
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, 0.2, 0.3, 0.4), center = 0),
    "penalty.*single number.*same length"
  )
})

test_that("check_penalty_shape rejects center with wrong length", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = 0.5, center = c(0, 0)),
    "center.*single number.*same length"
  )
  expect_error(
    check_penalty_shape(theta, penalty = 0.5, center = c(0, 0, 0, 0)),
    "center.*single number.*same length"
  )
})

test_that("check_penalty_shape rejects negative penalty values", {
  theta <- c(1, 2, 3)
  expect_error(
    check_penalty_shape(theta, penalty = -1, center = 0),
    "penalty.*non-negative"
  )
  expect_error(
    check_penalty_shape(theta, penalty = c(0.1, -0.2, 0.3), center = 0),
    "penalty.*non-negative"
  )
})

test_that("check_penalty_shape accepts zero penalty values", {
  theta <- c(1, 2, 3)
  expect_no_error(check_penalty_shape(theta, penalty = 0, center = 0))
  expect_no_error(check_penalty_shape(theta, penalty = c(0, 0, 0), center = 0))
})

# check_survival_data_valid -----------------------------------------------

test_that("check_survival_data_valid accepts valid survival data", {
  expect_no_error(check_survival_data_valid(delta = c(0, 1, 1, 0), time = c(1, 2, 3, 4)))
  expect_no_error(check_survival_data_valid(delta = c(1, 1, 1), time = c(0.5, 1.5, 10)))
  expect_no_error(check_survival_data_valid(delta = c(0, 0, 0), time = c(1, 2, 3)))
})

test_that("check_survival_data_valid accepts data with NAs", {
  # NAs in delta should be ignored

  expect_no_error(check_survival_data_valid(delta = c(0, 1, NA, 0), time = c(1, 2, 3, 4)))
  # NAs in time should be ignored
  expect_no_error(check_survival_data_valid(delta = c(0, 1, 1, 0), time = c(1, NA, 3, 4)))
  # NAs in both

  expect_no_error(check_survival_data_valid(delta = c(0, NA, 1, 0), time = c(1, NA, 3, 4)))
})

test_that("check_survival_data_valid rejects delta values other than 0 or 1", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 2), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
  expect_error(
    check_survival_data_valid(delta = c(-1, 0, 1), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
  expect_error(
    check_survival_data_valid(delta = c(0.5, 1, 0), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
})

test_that("check_survival_data_valid rejects delta with invalid values even when NAs present", {
  expect_error(
    check_survival_data_valid(delta = c(0, NA, 2), time = c(1, 2, 3)),
    "event indicator.*zero or one"
  )
})

test_that("check_survival_data_valid rejects non-positive times", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(0, 2, 3)),
    "times.*non-negative"
  )
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(-1, 2, 3)),
    "times.*non-negative"
  )
})

test_that("check_survival_data_valid rejects non-positive times even when NAs present", {
  expect_error(
    check_survival_data_valid(delta = c(0, 1, 1), time = c(NA, -1, 3)),
    "times.*non-negative"
  )
})

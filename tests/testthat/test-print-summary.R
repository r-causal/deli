# Tests for print() and summary() S7 methods for MEstimator and GMMEstimator

# The S-value column (S = -log2(P)) is highly sensitive: P is often tiny, so a
# last-ULP difference in the sandwich standard error, whether from a different
# BLAS reduction order or a different platform, moves S by enough to flip the
# fourth displayed decimal when the value sits on a rounding boundary. The
# named-parameter mtcars fit lands at S = 19.126949..., a hair below the 19.12695
# boundary between 19.1269 and 19.1270. Structural snapshots normalize that
# column so formatting is checked without pinning the boundary digit; the S-value
# is verified numerically with tolerance below and in test-inference.R.
stabilize_svalue <- function(lines) {
  sub("[0-9]+\\.[0-9]{4}$", "<S-value>", lines)
}

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

make_fitted_named <- function() {
  y <- mtcars$mpg
  X <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(intercept = 0, wt = 0, hp = 0)
  )
  estimate(m)
}

# print() snapshot tests -------------------------------------------------------

test_that("print() unfitted MEstimator", {
  y <- c(1, 2, 3)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_snapshot(print(m))
})

test_that("print() fitted single-parameter MEstimator", {
  m <- make_fitted_mean()
  expect_snapshot(print(m))
})

test_that("print() fitted multi-parameter MEstimator", {
  m <- make_fitted_mean_variance()
  expect_snapshot(print(m))
})

test_that("print() fitted named-parameter MEstimator", {
  m <- make_fitted_named()
  expect_snapshot(print(m))
})

test_that("print() unfitted GMMEstimator", {
  y <- c(1, 2, 3)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  g <- GMMEstimator(stacked_equations = psi, init = c(mean = 0))
  expect_snapshot(print(g))
})

test_that("print() fitted GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  g <- GMMEstimator(stacked_equations = psi, init = c(mean = 0))
  g <- estimate(g)
  expect_snapshot(print(g))
})

test_that("print() returns object invisibly", {
  m <- make_fitted_mean()
  out <- capture.output({
    result <- print(m)
  })
  expect_identical(result, m)
})

# summary() tests --------------------------------------------------------------

test_that("summary() errors on unfitted MEstimator", {
  y <- c(1, 2, 3)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(summary(m))
})

test_that("summary() returns an EstimatorSummary with expected properties", {
  m <- make_fitted_mean()
  s <- summary(m)
  expect_s3_class(s, "deli::EstimatorSummary")
  expect_no_error(s@theta)
  expect_no_error(s@se)
  expect_no_error(s@ci)
  expect_no_error(s@z)
  expect_no_error(s@p)
  expect_no_error(s@s)
  expect_no_error(s@n_obs)
  expect_no_error(s@n_params)
  expect_no_error(s@alpha)
  expect_no_error(s@estimator)
})

test_that("summary() values match individual inference functions", {
  m <- make_fitted_mean()
  s <- summary(m)
  expect_equal(s@theta, m@theta)
  expect_equal(s@se, sqrt(diag(m@variance)))
  expect_equal(s@ci, confidence_intervals(m))
  expect_equal(s@z, z_scores(m))
  expect_equal(s@p, p_values(m))
  expect_equal(s@s, s_values(m))
  expect_equal(s@n_obs, m@n_obs)
  expect_equal(s@n_params, m@n_params)
})

test_that("summary() with custom alpha changes CIs", {
  m <- make_fitted_mean()
  s90 <- summary(m, alpha = 0.10)
  s95 <- summary(m, alpha = 0.05)
  width_90 <- s90@ci[1, "upper"] - s90@ci[1, "lower"]
  width_95 <- s95@ci[1, "upper"] - s95@ci[1, "lower"]
  expect_true(width_90 < width_95)
})

# print(summary()) snapshot tests ----------------------------------------------

test_that("print(summary()) single-parameter output", {
  m <- make_fitted_mean()
  expect_snapshot(print(summary(m)), transform = stabilize_svalue)
})

test_that("print(summary()) multi-parameter output", {
  m <- make_fitted_mean_variance()
  expect_snapshot(print(summary(m)), transform = stabilize_svalue)
})

test_that("print(summary()) named-parameter output", {
  m <- make_fitted_named()
  s <- summary(m)
  expect_snapshot(print(s), transform = stabilize_svalue)
  # The boundary-sensitive S-value is checked numerically instead of by the
  # normalized snapshot above.
  expect_equal(
    unname(s@s),
    c(270.510082, 31.231015, 19.126949),
    tolerance = 1e-5
  )
})

test_that("print(summary()) GMMEstimator output", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  g <- GMMEstimator(stacked_equations = psi, init = c(mean = 0))
  g <- estimate(g)
  expect_snapshot(print(summary(g)), transform = stabilize_svalue)
})

test_that("print(summary()) with alpha = 0.10", {
  m <- make_fitted_mean()
  expect_snapshot(print(summary(m, alpha = 0.10)), transform = stabilize_svalue)
})

# summary() subset tests -------------------------------------------------------

test_that("summary() honors subset and keeps the selected parameter rows", {
  m <- make_fitted_named()
  full <- summary(m)
  sub <- summary(m, subset = c(1, 3))
  expect_equal(sub@theta, full@theta[c(1, 3)])
  expect_equal(sub@se, full@se[c(1, 3)])
  expect_equal(sub@ci, full@ci[c(1, 3), , drop = FALSE])
  expect_equal(sub@z, full@z[c(1, 3)])
  expect_equal(sub@p, full@p[c(1, 3)])
  expect_equal(sub@s, full@s[c(1, 3)])
  # n_params still reports the full model size, matching Python print_results
  expect_equal(sub@n_params, full@n_params)
})

test_that("summary() without subset returns every parameter", {
  m <- make_fitted_named()
  full <- summary(m)
  expect_length(full@theta, 3L)
  expect_equal(full@theta, m@theta)
})

test_that("summary() errors on an out-of-range subset", {
  m <- make_fitted_named()
  expect_error(summary(m, subset = c(1, 5)), "between 1 and 3")
})

test_that("summary() subset row labels track the original parameter index", {
  m <- make_fitted_mean_variance()
  s <- summary(m, subset = 2)
  expect_named(s@theta, "theta_2")
  expect_snapshot(print(s), transform = stabilize_svalue)
})

test_that("print(summary()) honors subset", {
  m <- make_fitted_named()
  expect_snapshot(
    print(summary(m, subset = c(1, 3))),
    transform = stabilize_svalue
  )
})

# print() estimator subset tests -----------------------------------------------

test_that("print() honors subset on the coefficient list", {
  m <- make_fitted_named()
  expect_snapshot(print(m, subset = c(1, 3)))
})

test_that("print() errors on an out-of-range subset", {
  m <- make_fitted_named()
  expect_error(print(m, subset = c(1, 5)), "between 1 and 3")
})

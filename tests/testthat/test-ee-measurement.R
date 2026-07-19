# Tests for measurement error estimating equations

# ee_rogan_gladen tests -------------------------------------------------------

test_that("ee_rogan_gladen returns 4-by-n matrix", {
  set.seed(42)
  n <- 100
  # External validation (r=0): have both y and y_star
  # Main study (r=1): only have y_star
  r <- c(rep(0, 40), rep(1, 60))
  y_true <- rbinom(n, 1, 0.4)
  # Imperfect measurement: sens=0.9, spec=0.8
  y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.8))
  # y is only observed in validation (r=0), fill rest with 0
  y <- ifelse(r == 0, y_true, 0)

  theta <- c(0.4, 0.4, 0.9, 0.8)
  result <- ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), n)
})

test_that("ee_rogan_gladen solves via MEstimator", {
  set.seed(123)
  n <- 500
  true_prev <- 0.3
  sens_true <- 0.9
  spec_true <- 0.85

  y_true <- rbinom(n, 1, true_prev)
  y_star <- ifelse(
    y_true == 1,
    rbinom(n, 1, sens_true),
    1 - rbinom(n, 1, spec_true)
  )
  r <- c(rep(0, 200), rep(1, 300))
  y <- ifelse(r == 0, y_true, 0)

  psi <- function(theta) {
    ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.3, 0.3, 0.8, 0.8))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Sensitivity and specificity should be in (0, 1)
  expect_true(m@theta[3] > 0 && m@theta[3] < 1)
  expect_true(m@theta[4] > 0 && m@theta[4] < 1)
})

test_that("ee_rogan_gladen with perfect measurement recovers naive mean", {
  set.seed(42)
  n <- 300
  y_true <- rbinom(n, 1, 0.5)
  y_star <- y_true # perfect measurement
  r <- c(rep(0, 100), rep(1, 200))
  y <- ifelse(r == 0, y_true, 0)

  psi <- function(theta) {
    ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.5, 0.5, 0.9, 0.9))
  m <- estimate(m, solver = "nleqslv")

  # With perfect measurement, corrected and naive should be close
  expect_equal(unname(m@theta[1]), unname(m@theta[2]), tolerance = 0.05)
  # Sensitivity and specificity should be ~1
  expect_equal(unname(m@theta[3]), 1, tolerance = 0.05)
  expect_equal(unname(m@theta[4]), 1, tolerance = 0.05)
})

# ee_rogan_gladen_extended tests -----------------------------------------------

test_that("ee_rogan_gladen_extended returns correct shape", {
  set.seed(42)
  n <- 100
  r <- c(rep(0, 40), rep(1, 60))
  y_true <- rbinom(n, 1, 0.4)
  y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.8))
  y <- ifelse(r == 0, y_true, 0)
  X <- cbind(rep(1, n)) # intercept only

  p <- ncol(X)
  theta <- c(0.4, rep(0, p), rep(0, p))
  result <- ee_rogan_gladen_extended(
    theta,
    y = y,
    y_star = y_star,
    r = r,
    X = X
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1 + 2 * p)
  expect_equal(ncol(result), n)
})

test_that("ee_rogan_gladen_extended solves via MEstimator", {
  set.seed(123)
  n <- 500
  y_true <- rbinom(n, 1, 0.3)
  y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
  r <- c(rep(0, 200), rep(1, 300))
  y <- ifelse(r == 0, y_true, 0)
  X <- cbind(rep(1, n))

  psi <- function(theta) {
    ee_rogan_gladen_extended(theta, y = y, y_star = y_star, r = r, X = X)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0.5, 1, 1))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  # Corrected proportion should be between 0 and 1
  expect_true(m@theta[1] > 0 && m@theta[1] < 1)
})

# ee_regression_calibration tests ---------------------------------------------

test_that("ee_regression_calibration returns correct shape", {
  set.seed(42)
  n <- 100
  r <- c(rep(0, 40), rep(1, 60))
  a_true <- rbinom(n, 1, 0.5)
  a_star <- ifelse(a_true == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.1))
  a <- ifelse(r == 0, a_true, 0)

  theta <- c(1, 0.5, 0.5) # corrected coef + 2 calibration params
  result <- ee_regression_calibration(
    theta,
    beta = 0.5,
    a = a,
    a_star = a_star,
    r = r
  )
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

test_that("ee_regression_calibration returns 2 + ncol(X) rows with a calibration X", {
  set.seed(43)
  n <- 100
  r <- c(rep(0, 40), rep(1, 60))
  a_true <- rbinom(n, 1, 0.5)
  a_star <- ifelse(a_true == 1, rbinom(n, 1, 0.9), rbinom(n, 1, 0.1))
  a <- ifelse(r == 0, a_true, 0)
  X <- cbind(1, rnorm(n))

  # corrected coefficient, then the calibration coefficients for
  # cbind(a_star, X): 1 + ncol(X) of them, so length 2 + ncol(X).
  theta <- rep(0, 2 + ncol(X))
  result <- ee_regression_calibration(
    theta,
    beta = 0.5,
    a = a,
    a_star = a_star,
    r = r,
    X = X
  )
  expect_equal(nrow(result), 2 + ncol(X))
  expect_equal(ncol(result), n)
})

test_that("ee_regression_calibration solves via MEstimator", {
  set.seed(789)
  n <- 500
  r <- c(rep(0, 200), rep(1, 300))
  a_true <- rbinom(n, 1, 0.5)
  a_star <- ifelse(a_true == 1, rbinom(n, 1, 0.85), rbinom(n, 1, 0.1))
  a <- ifelse(r == 0, a_true, 0)

  psi <- function(theta) {
    ee_regression_calibration(theta, beta = 0.8, a = a, a_star = a_star, r = r)
  }

  m <- MEstimator(stacked_equations = psi, init = c(1, 0.1, 0.5))
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
})

# NA masking parity ------------------------------------------------------------
#
# In validation designs the gold-standard measurement is only observed in the
# external sample (r == 0) and is legitimately missing in the main study
# (r == 1). Python masks that measurement with a placeholder before multiplying
# by the sample indicator, so NA values in the masked positions do not
# propagate. These tests pin the same tolerance for the R port. The fixtures are
# generated by
# tests/testthat/fixtures/generate_measurement_na_fixtures.py, which fits the
# Cole et al. (2023) examples with literal NaN in the masked positions. NA
# positions are serialized as JSON null and read back as NA.

# ee_rogan_gladen ----

test_that("ee_rogan_gladen returns finite rows with NA main-study y", {
  ref <- load_fixture("ee_rogan_gladen_na")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r

  # The gold-standard y is NA exactly in the main study (r == 1)
  expect_true(anyNA(y[r == 1]))
  expect_false(anyNA(y[r == 0]))

  ef <- ee_rogan_gladen(ref$theta, y = y, y_star = y_star, r = r)
  expect_true(all(is.finite(ef)))
})

test_that("ee_rogan_gladen fits and matches Python with NA main-study y", {
  ref <- load_fixture("ee_rogan_gladen_na")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r

  psi <- function(theta) {
    ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(m@variance)))
  expect_python_match(m, "ee_rogan_gladen_na", tolerance = 1e-5)
})

test_that("ee_rogan_gladen propagates NA outside the masked positions", {
  ref <- load_fixture("ee_rogan_gladen_na")
  y_star <- ref$y_star
  r <- ref$r

  # Clean main-study placeholder, then inject NA at a validation row (r == 0),
  # a position Python does not mask. The NA must still propagate so that the
  # eventual fix does not over-mask.
  y <- ref$y
  y[r == 1] <- 0
  y[which(r == 0)[1]] <- NA

  ef <- ee_rogan_gladen(ref$theta, y = y, y_star = y_star, r = r)
  expect_true(!all(is.finite(ef)))
})

# ee_rogan_gladen_extended ----

test_that("ee_rogan_gladen_extended returns finite rows with NA main-study y", {
  ref <- load_fixture("ee_rogan_gladen_extended_na")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r
  X <- ref$X

  expect_true(anyNA(y[r == 1]))
  expect_false(anyNA(y[r == 0]))

  ef <- ee_rogan_gladen_extended(
    ref$theta,
    y = y,
    y_star = y_star,
    r = r,
    X = X
  )
  expect_true(all(is.finite(ef)))
})

test_that("ee_rogan_gladen_extended fits and matches Python with NA main-study y", {
  ref <- load_fixture("ee_rogan_gladen_extended_na")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r
  X <- ref$X

  psi <- function(theta) {
    ee_rogan_gladen_extended(theta, y = y, y_star = y_star, r = r, X = X)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(m@variance)))
  expect_python_match(m, "ee_rogan_gladen_extended_na", tolerance = 1e-5)
})

test_that("ee_rogan_gladen_extended propagates NA outside the masked positions", {
  ref <- load_fixture("ee_rogan_gladen_extended_na")
  y_star <- ref$y_star
  r <- ref$r
  X <- ref$X

  y <- ref$y
  y[r == 1] <- 0
  y[which(r == 0)[1]] <- NA

  ef <- ee_rogan_gladen_extended(
    ref$theta,
    y = y,
    y_star = y_star,
    r = r,
    X = X
  )
  expect_true(!all(is.finite(ef)))
})

# ee_regression_calibration ----

test_that("ee_regression_calibration returns finite rows with NA main-study a", {
  ref <- load_fixture("ee_regression_calibration_na")
  a <- ref$a
  a_star <- ref$a_star
  r <- ref$r

  # The gold-standard a is NA exactly in the main study (r == 1)
  expect_true(anyNA(a[r == 1]))
  expect_false(anyNA(a[r == 0]))

  # theta = [corrected coef, 2 calibration params, 2 naive-model params]
  theta_calib <- ref$theta[1:3]
  beta <- ref$theta[5]
  ef <- ee_regression_calibration(
    theta_calib,
    beta = beta,
    a = a,
    a_star = a_star,
    r = r
  )
  expect_true(all(is.finite(ef)))
})

test_that("ee_regression_calibration fits and matches Python with NA main-study a", {
  ref <- load_fixture("ee_regression_calibration_na")
  a <- ref$a
  a_star <- ref$a_star
  y <- ref$y
  r <- ref$r

  # Stacked: calibration (3 params) + naive logistic (2 params)
  X_main <- cbind(1, a_star)

  psi <- function(theta) {
    theta_calib <- theta[1:3]
    theta_main <- theta[4:5]

    # Naive logistic (main study only); the naive outcome is masked in the
    # external data before use
    y_safe <- ifelse(r == 0, -999, y)
    ee_logit <- ee_regression(
      theta_main,
      X = X_main,
      y = y_safe,
      model = "logistic"
    )
    ee_logit <- ee_logit * rep(r, each = nrow(ee_logit))

    ee_calib <- ee_regression_calibration(
      theta_calib,
      beta = theta_main[2],
      a = a,
      a_star = a_star,
      r = r
    )

    rbind(ee_calib, ee_logit)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_true(all(is.finite(m@theta)))
  expect_true(all(is.finite(m@variance)))
  expect_python_match(m, "ee_regression_calibration_na", tolerance = 1e-5)
})

test_that("ee_regression_calibration propagates NA outside the masked positions", {
  ref <- load_fixture("ee_regression_calibration_na")
  a_star <- ref$a_star
  r <- ref$r

  # Clean main-study placeholder, then inject NA at a validation row (r == 0),
  # where the gold-standard a is required. The NA must still propagate.
  a <- ref$a
  a[r == 1] <- 0
  a[which(r == 0)[1]] <- NA

  theta_calib <- ref$theta[1:3]
  beta <- ref$theta[5]
  ef <- ee_regression_calibration(
    theta_calib,
    beta = beta,
    a = a,
    a_star = a_star,
    r = r
  )
  expect_true(!all(is.finite(ef)))
})

# Weighted extended Rogan-Gladen parity -----------------------------
#
# The extended correction fits logistic sensitivity and specificity nuisance
# models. Those models must receive the observation weights, not just the
# corrected-mean row. The fixture is generated by
# tests/testthat/fixtures/generate_rogan_gladen_extended_weighted_fixture.py,
# which fits the Cole et al. (2023) example with non-constant weights in Python
# Delicatessen.

test_that("ee_rogan_gladen_extended scales the nuisance rows by weights", {
  set.seed(7)
  n <- 12
  r <- c(rep(0, 6), rep(1, 6))
  y <- c(0, 1, 0, 1, 1, 0, rep(0, 6))
  y_star <- c(0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0)
  X <- cbind(rep(1, n)) # intercept only, so p = 1
  weights <- c(2, 1, 3, 1, 2, 1, 1, 1, 1, 2, 1, 3)
  theta <- c(0.5, 0.2, -0.1)

  ef_w <- ee_rogan_gladen_extended(
    theta,
    y = y,
    y_star = y_star,
    r = r,
    X = X,
    weights = weights
  )
  ef_u <- ee_rogan_gladen_extended(
    theta,
    y = y,
    y_star = y_star,
    r = r,
    X = X
  )

  # The sensitivity (row 2) and specificity (row 3) nuisance scores are linear
  # in the weights, so the weighted rows are the unweighted rows scaled by the
  # per-observation weight.
  expect_equal(ef_w[2, ], weights * ef_u[2, ])
  expect_equal(ef_w[3, ], weights * ef_u[3, ])
})

test_that("ee_rogan_gladen_extended fits and matches Python with weights", {
  ref <- load_fixture("ee_rogan_gladen_extended_weighted")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r
  X <- as.matrix(ref$X)
  weights <- ref$weights

  psi <- function(theta) {
    ee_rogan_gladen_extended(
      theta,
      y = y,
      y_star = y_star,
      r = r,
      X = X,
      weights = weights
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "lm")

  expect_python_match(m, "ee_rogan_gladen_extended_weighted", tolerance = 1e-5)
})

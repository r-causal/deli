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

# confidence_intervals() -------------------------------------------------------

test_that("confidence_intervals() returns a matrix with 2 columns", {
  m <- make_fitted_mean()
  ci <- confidence_intervals(m)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  expect_equal(nrow(ci), length(m@theta))
  expect_equal(colnames(ci), c("lower", "upper"))
})

test_that("confidence_intervals() matches Python reference values", {
  ref <- load_fixture("ee_mean")
  m <- make_fitted_mean()
  ci <- confidence_intervals(m, alpha = 0.05)
  expect_equal(unname(ci[, "lower"]), ref$ci_lower, tolerance = 1e-4)
  expect_equal(unname(ci[, "upper"]), ref$ci_upper, tolerance = 1e-4)
})

test_that("confidence_intervals() with alpha = 0.10 gives narrower CIs", {
  m <- make_fitted_mean()
  ci_default <- confidence_intervals(m, alpha = 0.05)
  ci_narrow <- confidence_intervals(m, alpha = 0.10)
  # Narrower CI means larger lower bound and smaller upper bound

  expect_gt(ci_narrow[1, "lower"], ci_default[1, "lower"])
  expect_lt(ci_narrow[1, "upper"], ci_default[1, "upper"])
})

test_that("confidence_intervals() errors before estimation", {
  y <- c(1, 2, 3)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(confidence_intervals(m))
})

# z_scores() -------------------------------------------------------------------

test_that("z_scores() returns numeric vector same length as theta", {
  m <- make_fitted_mean()
  z <- z_scores(m)
  expect_true(is.numeric(z))
  expect_length(z, length(m@theta))
})

test_that("z_scores() equals theta / SE when null = 0", {
  m <- make_fitted_mean()
  z <- z_scores(m, null = 0)
  se <- sqrt(diag(m@variance))
  expected_z <- m@theta / se
  expect_equal(z, expected_z)
})

test_that("z_scores() with custom null value works", {
  m <- make_fitted_mean()
  null_val <- 4
  z <- z_scores(m, null = null_val)
  se <- sqrt(diag(m@variance))
  expected_z <- (m@theta - null_val) / se
  expect_equal(z, expected_z)
})

# p_values() -------------------------------------------------------------------

test_that("p_values() returns numeric vector between 0 and 1", {
  m <- make_fitted_mean()
  p <- p_values(m)
  expect_true(is.numeric(p))
  expect_length(p, length(m@theta))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("p_values() matches 2 * pnorm(-abs(z))", {
  m <- make_fitted_mean()
  z <- z_scores(m, null = 0)
  p <- p_values(m, null = 0)
  expected_p <- 2 * pnorm(-abs(z))
  expect_equal(p, expected_p)
})

# s_values() -------------------------------------------------------------------

test_that("s_values() returns -log2(p_value)", {
  m <- make_fitted_mean()
  s <- s_values(m, null = 0)
  p <- p_values(m, null = 0)
  expected_s <- -log2(p)
  expect_equal(s, expected_s)
})

test_that("s_values() are always non-negative", {
  m <- make_fitted_mean()
  s <- s_values(m, null = 0)
  expect_true(all(s >= 0))
})

# An underflowed P-value -------------------------------------------------------
#
# The largest surprisal a double can name is 1074 bits, because the smallest
# positive double is 2^-1074. A P-value below that is not representable and
# arrives as exactly zero, so S = -log2(0) is Inf. That is the limit the
# surprisal is heading toward and the fit really does carry more bits than a
# double can hold, so the value is reported as Inf rather than capped at the
# representable maximum, which would name a number the evidence does not
# support. These tests pin that.

make_fitted_p_underflow <- function() {
  # One hundred observations alternating 9 and 11: the mean is 10 with a
  # sandwich standard error of 0.1, so Z is about 100 and the two-sided normal
  # tail probability, around 1e-2174, underflows to zero. The surprisal it
  # stands for is roughly 7220 bits against the 1074 a double can name.
  y <- rep(c(9, 11), 50)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  estimate(m)
}

test_that("s_values() reports Inf for a P-value that underflows to zero", {
  m <- make_fitted_p_underflow()
  # The Z-score itself is finite, so the zero is an underflowed tail
  # probability rather than a degenerate standard error.
  expect_equal(unname(z_scores(m)), 100, tolerance = 1e-6)
  expect_identical(unname(p_values(m)), 0)
  expect_identical(unname(s_values(m)), Inf)
})

test_that("s_values() stays finite for a P-value that does not underflow", {
  m <- make_fitted_mean()
  p <- p_values(m)
  expect_true(all(p > 0))
  expect_equal(s_values(m), -log2(p))
  expect_true(all(is.finite(s_values(m))))
})

test_that("tidy() reports an underflowed S-value as Inf", {
  m <- make_fitted_p_underflow()
  td <- tidy(m)
  expect_identical(td$p.value, 0)
  expect_identical(td$s.value, Inf)
})

# null length validation ------------------------------------------------------
#
# Python raises a broadcast ValueError when null has a length that is neither 1
# nor the parameter count. R instead recycles: a length-3 null on a 2-parameter
# fit warns twice about recycling but returns three z-scores, and those flow into
# p_values and s_values. These tests pin an informative abort for a mismatched
# null length across all three utilities, and confirm that a scalar null and a
# per-parameter null still work.

test_that("z_scores() aborts on a null length that is neither 1 nor the parameter count", {
  m <- make_fitted_mean_variance()
  expect_error(
    z_scores(m, null = c(0, 0, 0)),
    regexp = "null"
  )
})

test_that("p_values() aborts on a mismatched null length", {
  m <- make_fitted_mean_variance()
  expect_error(
    p_values(m, null = c(0, 0, 0)),
    regexp = "null"
  )
})

test_that("s_values() aborts on a mismatched null length", {
  m <- make_fitted_mean_variance()
  expect_error(
    s_values(m, null = c(0, 0, 0)),
    regexp = "null"
  )
})

test_that("z_scores() accepts a scalar null and a per-parameter null", {
  m <- make_fitted_mean_variance()
  se <- sqrt(diag(m@variance))
  expect_equal(z_scores(m, null = 0), m@theta / se)
  expect_equal(
    z_scores(m, null = c(1, 2)),
    (m@theta - c(1, 2)) / se
  )
})

# ---- deli_polygamma() -------------------------------------------------------

test_that("deli_polygamma() matches psigamma() for n=0 (digamma)", {
  x_vals <- c(0.5, 1, 2, 5, 10, 100)
  for (x in x_vals) {
    expect_equal(deli_polygamma(0, x), psigamma(x, deriv = 0))
  }
})

test_that("deli_polygamma() matches psigamma() for n=1 (trigamma)", {
  x_vals <- c(0.5, 1, 2, 5, 10)
  for (x in x_vals) {
    expect_equal(deli_polygamma(1, x), psigamma(x, deriv = 1))
  }
})

test_that("deli_polygamma() matches psigamma() for higher-order derivatives", {
  x_vals <- c(1, 2, 5)
  for (n in 2:4) {
    for (x in x_vals) {
      expect_equal(deli_polygamma(n, x), psigamma(x, deriv = n))
    }
  }
})

test_that("deli_polygamma() works with vector input", {
  x_vec <- c(0.5, 1, 2, 5, 10)
  expect_equal(deli_polygamma(0, x_vec), psigamma(x_vec, deriv = 0))
  expect_equal(deli_polygamma(1, x_vec), psigamma(x_vec, deriv = 1))
  expect_equal(deli_polygamma(2, x_vec), psigamma(x_vec, deriv = 2))
})

test_that("deli_polygamma() handles large x values", {
  expect_equal(deli_polygamma(0, 1000), psigamma(1000, deriv = 0))
  expect_equal(deli_polygamma(1, 1000), psigamma(1000, deriv = 1))
})

# deli_polygamma(n, x) takes the order first, following Python's
# scipy.special.polygamma; psigamma(x, deriv = n) takes it second. A positional
# find-and-replace between the two silently computes a different quantity, so the
# correspondence is pinned in both directions with two distinct arguments, and
# the swapped result is asserted to disagree.

test_that("deli_polygamma() takes the order first where psigamma() takes it second", {
  expect_equal(deli_polygamma(1, 3), psigamma(3, deriv = 1))
  expect_equal(deli_polygamma(3, 1), psigamma(1, deriv = 3))
  expect_false(isTRUE(all.equal(
    psigamma(3, deriv = 1),
    psigamma(1, deriv = 3)
  )))
})

# ---- deli_digamma() ---------------------------------------------------------

test_that("deli_digamma() matches base::digamma() for positive values", {
  z_vals <- c(0.5, 1, 2, 5, 10, 100)
  for (z in z_vals) {
    expect_equal(deli_digamma(z), digamma(z))
  }
})

test_that("deli_digamma() matches base::digamma() for negative non-integer values", {
  z_vals <- c(-0.5, -1.5, -2.7)
  for (z in z_vals) {
    expect_equal(deli_digamma(z), digamma(z))
  }
})

test_that("deli_digamma() works with vector input", {
  z_vec <- c(0.5, 1, 2, 5, 10)
  expect_equal(deli_digamma(z_vec), digamma(z_vec))
})

test_that("deli_digamma() handles large values", {
  expect_equal(deli_digamma(1000), digamma(1000))
  expect_equal(deli_digamma(1e6), digamma(1e6))
})

test_that("deli_digamma() at negative integers returns NaN without warning", {
  # digamma is undefined at non-positive integers; deli_digamma handles
  # this gracefully by returning NaN without triggering a warning
  result0 <- expect_no_warning(deli_digamma(0))
  expect_true(is.nan(result0))
  result_neg1 <- expect_no_warning(deli_digamma(-1))
  expect_true(is.nan(result_neg1))
})

# ---- NA / NaN propagation ----------------------------------------
#
# The pole test `z <= 0 & z == trunc(z)` yields NA for a missing input, which
# poisons both the if() guard and the subscripted assignment: deli_digamma(NA)
# and deli_digamma(NaN) crash, and a vector with a NaN entry crashes with "NAs
# are not allowed in subscripted assignments". base::digamma and scipy both
# propagate NA/NaN, so deli_digamma must too, preserving the NA vs NaN
# distinction and keeping the documented NaN at the poles.

test_that("deli_digamma() propagates a scalar NA or NaN without error", {
  r_na <- expect_no_error(deli_digamma(NA_real_))
  expect_true(is.na(r_na))
  expect_false(is.nan(r_na))

  r_nan <- expect_no_error(deli_digamma(NaN))
  expect_true(is.nan(r_nan))
})

test_that("deli_digamma() propagates NA and NaN elementwise in a vector", {
  r <- expect_no_error(deli_digamma(c(1, NaN, 2)))
  expect_equal(r[1], digamma(1))
  expect_true(is.nan(r[2]))
  expect_equal(r[3], digamma(2))

  r2 <- deli_digamma(c(1, NA, 2))
  expect_equal(r2[1], digamma(1))
  expect_true(is.na(r2[2]))
  expect_false(is.nan(r2[2]))
  expect_equal(r2[3], digamma(2))
})

test_that("deli_digamma() keeps NaN at the poles alongside missing input", {
  r <- deli_digamma(c(NaN, 0, -2, 2.5))
  expect_true(is.nan(r[1])) # NaN input propagates
  expect_true(is.nan(r[2])) # pole at 0
  expect_true(is.nan(r[3])) # pole at -2
  expect_equal(r[4], digamma(2.5)) # regular value
})

test_that("deli_digamma() differs from base::digamma() only in the pole warning", {
  # The reference documentation states that the suppressed pole warning is the
  # only difference between the two, so every returned value has to agree: the
  # regular values, the NaN at each pole, and the NA and NaN propagated from
  # missing input. The comparison is `expect_identical()` rather than
  # `expect_equal()` because edition 3 gives `expect_equal()` a default numeric
  # tolerance, and waldo's tolerant path treats NA_real_ and NaN as equal, which
  # would leave the missing-value axis unchecked. The expected base R warning is
  # caught explicitly rather than suppressed.
  z <- c(-2.7, -1.5, 0.5, 1, 2.5, 10, 0, -1, -3, NA, NaN)
  expect_warning(base_values <- digamma(z), "NaNs produced")
  deli_values <- expect_no_warning(deli_digamma(z))
  expect_identical(deli_values, base_values)
})

# ---- standard_normal_cdf() --------------------------------------------------

test_that("standard_normal_cdf() matches pnorm() for various values", {
  x_vals <- c(-3, -2, -1, 0, 1, 2, 3)
  for (x in x_vals) {
    expect_equal(standard_normal_cdf(x), pnorm(x))
  }
})

test_that("standard_normal_cdf(0) returns 0.5", {
  expect_equal(standard_normal_cdf(0), 0.5)
})

test_that("standard_normal_cdf() works with vector input", {
  x_vec <- c(-3, -2, -1, 0, 1, 2, 3)
  expect_equal(standard_normal_cdf(x_vec), pnorm(x_vec))
})

test_that("standard_normal_cdf() handles extreme values", {
  # Very large positive -> 1
  expect_equal(standard_normal_cdf(10), pnorm(10))
  expect_equal(standard_normal_cdf(38), pnorm(38))
  # Very large negative -> 0
  expect_equal(standard_normal_cdf(-10), pnorm(-10))
  expect_equal(standard_normal_cdf(-38), pnorm(-38))
})

test_that("standard_normal_cdf() is monotonically increasing", {
  x_vec <- seq(-5, 5, by = 0.5)
  result <- standard_normal_cdf(x_vec)
  expect_true(all(diff(result) > 0))
})

test_that("standard_normal_cdf() is symmetric around 0", {
  x_vals <- c(0.5, 1, 2, 3)
  for (x in x_vals) {
    expect_equal(
      standard_normal_cdf(x) + standard_normal_cdf(-x),
      1
    )
  }
})

# ---- standard_normal_pdf() --------------------------------------------------

test_that("standard_normal_pdf() matches dnorm() for various values", {
  x_vals <- c(-3, -2, -1, 0, 1, 2, 3)
  for (x in x_vals) {
    expect_equal(standard_normal_pdf(x), dnorm(x))
  }
})

test_that("standard_normal_pdf(0) returns the mode density", {
  expect_equal(standard_normal_pdf(0), 1 / sqrt(2 * pi))
})

test_that("standard_normal_pdf() works with vector input", {
  x_vec <- c(-3, -2, -1, 0, 1, 2, 3)
  expect_equal(standard_normal_pdf(x_vec), dnorm(x_vec))
})

test_that("standard_normal_pdf() handles extreme values", {
  expect_equal(standard_normal_pdf(10), dnorm(10))
  expect_equal(standard_normal_pdf(-10), dnorm(-10))
  expect_equal(standard_normal_pdf(38), dnorm(38))
})

test_that("standard_normal_pdf() is symmetric around 0", {
  x_vals <- c(0.5, 1, 2, 3)
  for (x in x_vals) {
    expect_equal(standard_normal_pdf(x), standard_normal_pdf(-x))
  }
})

test_that("standard_normal_pdf() values are non-negative", {
  x_vec <- seq(-5, 5, by = 0.5)
  result <- standard_normal_pdf(x_vec)
  expect_true(all(result >= 0))
})

# ---- whole-vector inputs under exact autodiff --------------------
#
# Passing the whole parameter vector (a PrimalTangentVector under exact mode) to
# deli_polygamma, deli_digamma, standard_normal_cdf, and standard_normal_pdf must
# differentiate elementwise. Only is_pt / is_pt_array inputs were handled, so a
# PrimalTangentVector fell through to the numeric kernel and errored with
# "Non-numeric argument to mathematical function". The Jacobian is diagonal
# because each output depends on a single parameter, which discriminates a
# correct per-element routing from a mix-up. Cross-checks follow the autodiff
# convention: analytic at 1e-8, capprox at 1e-5.

test_that("standard_normal_cdf() differentiates a whole parameter vector under exact mode", {
  theta <- c(0.5, 1)
  f <- function(x) standard_normal_cdf(x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- diag(dnorm(theta)) # d/dx pnorm(x) = dnorm(x)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("standard_normal_pdf() differentiates a whole parameter vector under exact mode", {
  theta <- c(0.5, 1)
  f <- function(x) standard_normal_pdf(x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- diag(-theta * dnorm(theta)) # d/dx dnorm(x) = -x dnorm(x)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("deli_digamma() differentiates a whole parameter vector under exact mode", {
  theta <- c(2.5, 3.1)
  f <- function(x) deli_digamma(x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- diag(trigamma(theta)) # d/dx digamma(x) = trigamma(x)
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("deli_polygamma() differentiates a whole parameter vector under exact mode", {
  theta <- c(2.5, 3.1)
  f <- function(x) deli_polygamma(1, x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # d/dx polygamma(1, x) = polygamma(2, x) = psigamma(x, deriv = 2)
  analytic <- diag(psigamma(theta, deriv = 2))
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

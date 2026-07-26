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

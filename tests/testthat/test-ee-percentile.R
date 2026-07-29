# Tests for ee_percentile() and ee_positive_mean_deviation()

test_that("ee_percentile with q=0.5 returns correct shape", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  result <- suppressWarnings(ee_percentile(theta = c(5), y = y, q = 0.5))

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), length(y))
})

test_that("ee_percentile errors for invalid q", {
  y <- c(1, 2, 3, 4, 5)

  expect_error(
    ee_percentile(theta = c(3), y = y, q = 0),
    "must be between 0 and 1"
  )

  expect_error(
    ee_percentile(theta = c(3), y = y, q = 1),
    "must be between 0 and 1"
  )

  expect_error(
    ee_percentile(theta = c(3), y = y, q = -0.5),
    "must be between 0 and 1"
  )

  expect_error(
    ee_percentile(theta = c(3), y = y, q = 1.5),
    "must be between 0 and 1"
  )
})

test_that("ee_positive_mean_deviation returns 2-by-n matrix", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  result <- suppressWarnings(ee_positive_mean_deviation(theta = c(3, 5), y = y))

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), length(y))
})

# The percentile estimating function is non-smooth: its derivative is undefined
# at the solution, so the bread matrix does not exist for finite samples and the
# sandwich variance is not valid. ee_percentile() and ee_positive_mean_deviation()
# signal a non-differentiability warning on every direct call, matching Python
# delicatessen. The warning is advisory and does not change the returned
# contributions.
#
# Intended message for ee_percentile():
#   "The estimating equation is not differentiable at {.arg theta}. Therefore,
#    the bread matrix is not defined for finite samples, and the sandwich should
#    not be used to estimate the variance."
# Intended message for ee_positive_mean_deviation():
#   "The estimating equation for the median is not differentiable. Therefore,
#    the bread matrix is not defined for finite samples, and the sandwich should
#    not be used to estimate the variance."
# Both tests match on the stable substring "not differentiable".

test_that("ee_percentile warns that it is not differentiable", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  expect_warning(
    ee_percentile(theta = c(5), y = y, q = 0.5),
    "not differentiable"
  )
})

test_that("ee_percentile returns the correct contributions despite the warning", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  result <- suppressWarnings(ee_percentile(theta = c(5), y = y, q = 0.5))

  expected <- matrix(
    c(-0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5),
    nrow = 1
  )
  expect_equal(result, expected)
})

test_that("ee_positive_mean_deviation warns that the median is not differentiable", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  expect_warning(
    ee_positive_mean_deviation(theta = c(3, 5), y = y),
    "not differentiable"
  )
})

test_that("ee_positive_mean_deviation returns the correct contributions despite the warning", {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)

  result <- suppressWarnings(ee_positive_mean_deviation(theta = c(3, 5), y = y))

  expected <- rbind(
    c(-3, -3, -3, -3, -3, -1, 1, 3, 5),
    c(-0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5)
  )
  rownames(expected) <- c("positive_mean_deviation", "median")
  expect_equal(result, expected)
})

# The median row of ee_positive_mean_deviation() is a step function, so its
# Jacobian row is zero and no solver can search for the median. The documented
# advice is therefore to start theta[2] at the sample median and theta[1] at the
# positive mean deviation that matches it, under solver = "lm"; the example on
# that help page returns the sample median only because that start holds. Note
# that what is pinned below is minpack.lm's behavior, not deli's. A failure here
# means the dependency changed how it treats a singular Jacobian row, not that
# deli regressed. Even and odd n are both covered because the median row is
# exactly satisfiable only at even n without ties straddling the median.
test_that("the documented starting values hold the sample median under the lm solver", {
  skip_if_not_installed("minpack.lm")

  fit_at_documented_init <- function(y) {
    psi <- function(theta) {
      suppressWarnings(ee_positive_mean_deviation(theta, y = y))
    }
    init <- c(mean(2 * (y - median(y)) * (y > median(y))), median(y))
    m <- m_estimate(stacked_equations = psi, init = init, solver = "lm")
    unname(coef(m))
  }

  expect_equal(fit_at_documented_init(c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)), c(2, 3))

  set.seed(1)
  y_even <- stats::rexp(40)
  expect_equal(fit_at_documented_init(y_even)[2], median(y_even))

  set.seed(1)
  y_odd <- stats::rexp(41)
  expect_equal(fit_at_documented_init(y_odd)[2], median(y_odd))
})

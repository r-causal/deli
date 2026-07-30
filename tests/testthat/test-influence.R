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

# Parameter names --------------------------------------------------------------
#
# The columns of the returned matrix are the parameters, so they carry the
# parameter names as every other post-estimation accessor does. Most of the time
# they arrive on their own: estimate() names the bread and solve() propagates
# those names through the inverse. MASS::ginv() does not, so a singular bread is
# the case that needs the assignment, and the case these tests cover.

test_that("influence_functions() names its columns for the parameters", {
  m <- make_fitted_mean_variance()
  expect_equal(colnames(influence_functions(m)), names(coef(m)))
})

test_that("influence_functions() names its columns when the bread is singular", {
  skip_if_not_installed("MASS")
  # Duplicating a predictor makes the three-parameter bread rank two, so
  # solve() fails and the pseudo-inverse takes over. Solving only the intercept
  # and the wt slope leaves the duplicate frozen at zero, which keeps the
  # root-finder quiet while the bread stays rank deficient.
  d <- mtcars[c("mpg", "wt")]
  d$wt_copy <- d$wt
  X <- stats::model.matrix(mpg ~ wt + wt_copy, data = d)
  psi <- function(theta) {
    ee_regression(theta, X = X, y = d$mpg, model = "linear")
  }
  m <- m_estimate(
    stacked_equations = psi,
    init = c(`(Intercept)` = 0, wt = 0, wt_copy = 0),
    subset = c(1L, 2L)
  )
  # The fixture is only doing its job if solve() really cannot invert the bread.
  expect_error(solve(m@bread), "singular")

  expect_equal(
    colnames(influence_functions(m)),
    c("(Intercept)", "wt", "wt_copy")
  )
})

# ---- a clustered fit reports one row per group -------------------------------
#
# aggregate_efuncs() collapses the per-observation contributions within a group
# and labels the resulting columns with the group values. Those labels ride
# through the bread inverse and out of the transpose, so the rows of a clustered
# fit's influence functions are the groups rather than the observations, and they
# say which group each one is. Nothing else in the returned matrix says so: the
# row count is the number of groups and reads as a number of observations.

clustered_mean_fit <- function() {
  set.seed(42)
  n <- 40
  group <- rep(c("a", "b", "c", "d", "e"), each = 8)
  y <- stats::rnorm(5, sd = 2)[match(group, sort(unique(group)))] +
    stats::rnorm(n)
  psi <- function(theta) {
    aggregate_efuncs(ee_mean(theta, y = y), group = group)
  }
  list(fit = m_estimate(stacked_equations = psi, init = mean(y)), n = n)
}

test_that("influence_functions() labels the rows of a clustered fit by group", {
  case <- clustered_mean_fit()
  inf <- influence_functions(case$fit)

  expect_equal(nrow(inf), 5L)
  expect_lt(nrow(inf), case$n)
  expect_equal(rownames(inf), c("a", "b", "c", "d", "e"))
})

test_that("an unaggregated fit leaves the influence-function rows unlabeled", {
  set.seed(42)
  y <- stats::rnorm(20)
  fit <- m_estimate(
    stacked_equations = function(theta) ee_mean(theta, y = y),
    init = mean(y)
  )

  expect_equal(nrow(influence_functions(fit)), 20L)
  expect_null(rownames(influence_functions(fit)))
})

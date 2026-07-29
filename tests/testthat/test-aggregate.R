# Tests for aggregate_efuncs

test_that("aggregate_efuncs returns correct shape", {
  # 2 parameters, 6 observations, 3 groups
  ef <- matrix(1:12, nrow = 2, ncol = 6)
  group <- c(1, 1, 2, 2, 3, 3)

  result <- aggregate_efuncs(ef, group)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})

test_that("aggregate_efuncs sums within groups correctly", {
  # Simple case: 1 parameter, 4 observations, 2 groups
  ef <- matrix(c(1, 2, 3, 4), nrow = 1)
  group <- c("a", "a", "b", "b")

  result <- aggregate_efuncs(ef, group)
  expect_equal(result[1, 1], 3) # 1 + 2
  expect_equal(result[1, 2], 7) # 3 + 4
})

test_that("aggregate_efuncs with singletons returns same values", {
  ef <- matrix(c(10, 20, 30), nrow = 1)
  group <- c(1, 2, 3)

  result <- aggregate_efuncs(ef, group)
  expect_equal(ncol(result), 3)
  expect_equal(as.numeric(result), c(10, 20, 30))
})

test_that("aggregate_efuncs works with multiple parameters", {
  # 3 parameters, 4 observations, 2 groups
  ef <- matrix(
    c(
      1,
      2,
      3,
      4, # param 1
      5,
      6,
      7,
      8, # param 2
      9,
      10,
      11,
      12 # param 3
    ),
    nrow = 3,
    byrow = TRUE
  )
  group <- c(1, 1, 2, 2)

  result <- aggregate_efuncs(ef, group)
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), 2)
  expect_equal(result[1, 1], 3) # 1+2
  expect_equal(result[1, 2], 7) # 3+4
  expect_equal(result[2, 1], 11) # 5+6
  expect_equal(result[2, 2], 15) # 7+8
  expect_equal(result[3, 1], 19) # 9+10
  expect_equal(result[3, 2], 23) # 11+12
})

test_that("aggregate_efuncs treats a 1-D vector as one parameter across observations", {
  # Python (delicatessen.utilities.aggregate_efuncs) treats a length-n 1-D input
  # as a single parameter observed n times (shape 1-by-m), not n parameters
  # observed once. Pin that orientation so `as.matrix` on a vector cannot
  # transpose the intent into an n-by-1 column.
  ef <- c(1, 2, 3, 4)
  group <- c("a", "a", "b", "b")

  result <- aggregate_efuncs(ef, group)

  # Shape: one parameter, two groups (matches Python r.shape == (1, 2))
  expect_equal(dim(result), c(1L, 2L))
  # Values: 1 + 2 = 3 for group "a", 3 + 4 = 7 for group "b"
  expect_equal(result[1, 1], 3)
  expect_equal(result[1, 2], 7)
})

test_that("aggregate_efuncs 1-D input matches its explicit 1-row matrix form", {
  # A length-n vector and a 1-by-n matrix carry the same contract: one
  # parameter across n observations. Both must aggregate identically.
  values <- c(10, -1, 4, 5, -2, 3)
  group <- c(1, 1, 2, 3, 3, 3)

  from_vector <- aggregate_efuncs(values, group)
  from_matrix <- aggregate_efuncs(matrix(values, nrow = 1), group)

  expect_equal(dim(from_vector), c(1L, 3L))
  expect_equal(from_vector, from_matrix)
  # Group 1: 10 + -1 = 9; group 2: 4; group 3: 5 + -2 + 3 = 6
  expect_equal(as.numeric(from_vector), c(9, 4, 6))
})

test_that("aggregate_efuncs returns groups in sorted order", {
  # Python (delicatessen) uses np.unique on the group vector, which sorts the
  # unique groups. Columns must come back in sorted unique-group order
  # regardless of the order groups first appear in the input.
  ef <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 1)
  group <- c(2, 2, 3, 3, 1, 1)

  result <- aggregate_efuncs(ef, group)

  # Columns follow sorted groups 1, 2, 3, not first-appearance order 2, 3, 1.
  # Group 1: 5 + 6 = 11; group 2: 1 + 2 = 3; group 3: 3 + 4 = 7.
  expect_equal(dim(result), c(1L, 3L))
  expect_equal(as.numeric(result), c(11, 3, 7))
})

test_that("aggregate_efuncs sorts groups with multiple parameters", {
  # Two parameters observed across six units, groups out of order. Each
  # parameter's columns must be reordered to sorted group order together.
  ef <- matrix(
    c(
      1,
      2,
      3,
      4,
      5,
      6, # param 1
      10,
      20,
      30,
      40,
      50,
      60 # param 2
    ),
    nrow = 2,
    byrow = TRUE
  )
  group <- c(3, 3, 1, 1, 2, 2)

  result <- aggregate_efuncs(ef, group)

  # Sorted groups 1, 2, 3.
  # Group 1: cols 3,4 -> param1 3+4=7, param2 30+40=70
  # Group 2: cols 5,6 -> param1 5+6=11, param2 50+60=110
  # Group 3: cols 1,2 -> param1 1+2=3, param2 10+20=30
  expect_equal(dim(result), c(2L, 3L))
  expect_equal(result[1, ], c(7, 11, 3))
  expect_equal(result[2, ], c(70, 110, 30))
})

test_that("aggregate_efuncs errors on dimension mismatch", {
  ef <- matrix(1:6, nrow = 2)
  group <- c(1, 2) # length 2, but ef has 3 columns

  expect_error(aggregate_efuncs(ef, group), "group")
})

test_that("aggregate_efuncs with MEstimator gives cluster-robust variance", {
  set.seed(42)
  n <- 200
  n_groups <- 50
  group <- rep(1:n_groups, each = n / n_groups)

  # Generate clustered data: group-level intercept + noise
  group_effect <- rnorm(n_groups, sd = 2)
  y <- group_effect[group] + rnorm(n)

  # Without aggregation: standard variance
  psi_ind <- function(theta) {
    ee_mean(theta, y = y)
  }
  m_ind <- MEstimator(stacked_equations = psi_ind, init = mean(y))
  m_ind <- estimate(m_ind)

  # With aggregation: cluster-robust variance
  psi_agg <- function(theta) {
    ef <- ee_mean(theta, y = y)
    aggregate_efuncs(ef, group = group)
  }
  m_agg <- MEstimator(stacked_equations = psi_agg, init = mean(y))
  m_agg <- estimate(m_agg)

  # Point estimates should be the same
  expect_equal(m_ind@theta, m_agg@theta, tolerance = 1e-6)

  # Cluster-robust variance should be larger (clustered data)
  expect_true(m_agg@variance[1, 1] > m_ind@variance[1, 1])
})

test_that("aggregate_efuncs with unequal cluster sizes", {
  set.seed(123)
  # Unequal cluster sizes
  group <- c(rep(1, 10), rep(2, 5), rep(3, 15), rep(4, 20))
  n <- length(group)
  y <- rnorm(n, mean = 5)

  psi <- function(theta) {
    ef <- ee_mean(theta, y = y)
    aggregate_efuncs(ef, group = group)
  }

  m <- MEstimator(stacked_equations = psi, init = mean(y))
  m <- estimate(m)

  expect_true(is.finite(m@theta[1]))
  expect_true(m@variance[1, 1] > 0)
})

test_that("aggregate_efuncs works with regression EEs", {
  set.seed(789)
  n <- 100
  n_groups <- 25
  group <- rep(1:n_groups, each = n / n_groups)
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  X <- cbind(1, x)

  psi <- function(theta) {
    ef <- ee_regression(theta, X = X, y = y, model = "linear")
    aggregate_efuncs(ef, group = group)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 0))
  m <- estimate(m)

  expect_true(all(is.finite(m@theta)))
  expect_equal(unname(m@theta[1]), 1, tolerance = 0.5)
  expect_true(all(diag(m@variance) > 0))
})

# Parameter names --------------------------------------------------------------
#
# The rows of the aggregated matrix are the same parameters as the rows of the
# input, so labels on the input rows survive aggregation and reach the fitted
# estimator through the estimating function. The columns are a different matter:
# rowsum() labels its rows with the compact group indices rather than the
# original group values, so those labels are dropped instead of being passed off
# as group names.

test_that("aggregate_efuncs() keeps the row names of its input", {
  ef <- matrix(1:12, nrow = 2, dimnames = list(c("mu", "sigma2"), NULL))
  group <- c(1, 1, 2, 2, 3, 3)

  result <- aggregate_efuncs(ef, group)

  expect_equal(rownames(result), c("mu", "sigma2"))
})

test_that("aggregate_efuncs() leaves the aggregated columns unlabeled", {
  # The group values are 10, 20 and 30, so an integer index passed off as a
  # group label would be visibly wrong.
  ef <- matrix(1:12, nrow = 2, dimnames = list(c("mu", "sigma2"), NULL))
  group <- c(10, 10, 20, 20, 30, 30)

  result <- aggregate_efuncs(ef, group)

  expect_null(colnames(result))
})

test_that("aggregate_efuncs() returns an unnamed matrix for an unnamed input", {
  ef <- matrix(1:12, nrow = 2)
  group <- c(1, 1, 2, 2, 3, 3)

  expect_null(dimnames(aggregate_efuncs(ef, group)))
})

test_that("aggregate_efuncs() row names leave a fit undisturbed", {
  # The estimator labels its own matrices from `init`, so row names arriving on
  # the estimating function must reach the sandwich without displacing them.
  set.seed(4242)
  n <- 100
  group <- rep(1:25, each = 4)
  y <- rnorm(n, mean = 5)

  psi <- function(theta) {
    ef <- ee_mean(theta, y = y)
    rownames(ef) <- "labelled"
    aggregate_efuncs(ef, group = group)
  }

  m <- m_estimate(stacked_equations = psi, init = c(mu = 0))

  expect_equal(names(coef(m)), "mu")
  expect_equal(dimnames(m@meat), list("mu", "mu"))
  expect_true(m@variance[1, 1] > 0)
})

# ---- compute_meat() ----------------------------------------------------------

test_that("compute_meat() returns correct p-by-p matrix for known input", {
  # 2 parameters, 3 observations
  evaluations <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    ncol = 3
  )
  # Expected: tcrossprod(evaluations) = evaluations %*% t(evaluations)
  #   [1*1+3*3+5*5, 1*2+3*4+5*6]   [35, 44]
  #   [2*1+4*3+6*5, 2*2+4*4+6*6]   [44, 56]
  expected <- matrix(c(35, 44, 44, 56), nrow = 2, ncol = 2)
  result <- compute_meat(evaluations)
  expect_equal(dim(result), c(2L, 2L))
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("compute_meat() handles single observation (n=1)", {
  evaluations <- matrix(c(3, 5), nrow = 2, ncol = 1)
  expected <- matrix(c(9, 15, 15, 25), nrow = 2, ncol = 2)
  result <- compute_meat(evaluations)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("compute_meat() returns symmetric matrix", {
  set.seed(42)
  evaluations <- matrix(rnorm(20), nrow = 4, ncol = 5)
  result <- compute_meat(evaluations)
  expect_equal(result, t(result), tolerance = 1e-6)
})

test_that("compute_meat() returns 1x1 matrix for single parameter", {
  evaluations <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  # tcrossprod: sum(c(1,2,3)^2) = 14
  expected <- matrix(14, nrow = 1, ncol = 1)
  result <- compute_meat(evaluations)
  expect_equal(result, expected, tolerance = 1e-6)
})

# ---- build_sandwich() -------------------------------------------------------

test_that("build_sandwich() returns correct result for known matrices", {
  bread <- matrix(c(2, 0, 0, 3), nrow = 2, ncol = 2)
  meat <- matrix(c(4, 1, 1, 9), nrow = 2, ncol = 2)
  bread_inv <- solve(bread)
  expected <- bread_inv %*% meat %*% t(bread_inv)
  result <- build_sandwich(bread, meat)
  expect_equal(dim(result), c(2L, 2L))
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("build_sandwich() with identity bread returns meat", {
  meat <- matrix(c(4, 1, 1, 9), nrow = 2, ncol = 2)
  bread <- diag(2)
  result <- build_sandwich(bread, meat)
  expect_equal(result, meat, tolerance = 1e-6)
})

test_that("build_sandwich() with diagonal bread gives correct scaling", {
  bread <- diag(c(2, 4))
  meat <- diag(c(8, 16))
  expected <- diag(c(2, 1))
  result <- build_sandwich(bread, meat)
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("build_sandwich() with singular bread uses pseudo-inverse", {
  bread <- matrix(c(1, 2, 2, 4), nrow = 2, ncol = 2)
  meat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  result <- build_sandwich(bread, meat, allow_pinv = TRUE)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 2L))
})

test_that("build_sandwich() returns NULL when bread has NA", {
  bread <- matrix(c(1, NA, NA, 1), nrow = 2, ncol = 2)
  meat <- diag(2)
  result <- build_sandwich(bread, meat)
  expect_null(result)
})

test_that("build_sandwich() returns symmetric matrix for symmetric inputs", {
  bread <- matrix(c(3, 1, 1, 2), nrow = 2, ncol = 2)
  meat <- matrix(c(5, 2, 2, 4), nrow = 2, ncol = 2)
  result <- build_sandwich(bread, meat)
  expect_equal(result, t(result), tolerance = 1e-6)
})

test_that("build_sandwich() works for 1x1 matrices", {
  bread <- matrix(4, nrow = 1, ncol = 1)
  meat <- matrix(16, nrow = 1, ncol = 1)
  expected <- matrix(1, nrow = 1, ncol = 1)
  result <- build_sandwich(bread, meat)
  expect_equal(result, expected, tolerance = 1e-6)
})

# ---- compute_bread() --------------------------------------------------------

test_that("compute_bread() returns negative Jacobian for linear EE", {
  n <- 5
  func <- function(theta) {
    ee1 <- rep(1 - theta[1], n)
    ee2 <- rep(2 - theta[2], n)
    rbind(ee1, ee2)
  }
  theta <- c(1, 2)
  result <- compute_bread(func, theta)
  # Jacobian of sum(c_i - theta) = -n*I, bread = -(-n*I) = n*I
  expect_equal(dim(result), c(2L, 2L))
  expect_equal(unname(result), n * diag(2), tolerance = 1e-5)
})

test_that("compute_bread() returns correct dimensions", {
  n <- 10
  func <- function(theta) {
    matrix(rnorm(length(theta) * n), nrow = length(theta), ncol = n)
  }
  theta <- c(1, 2, 3)
  result <- compute_bread(func, theta)
  expect_equal(dim(result), c(3L, 3L))
})

# ---- finite_sample_correction() ---------------------------------------------

test_that("finite_sample_correction() with NULL returns meat unchanged", {
  meat <- matrix(c(4, 1, 1, 9), nrow = 2, ncol = 2)
  result <- finite_sample_correction(meat, n = 10, p = 2, adjustment = NULL)
  expect_equal(result, meat)
})

test_that("finite_sample_correction() applies HC1 correctly", {
  meat <- matrix(c(4, 1, 1, 9), nrow = 2, ncol = 2)
  n <- 10
  p <- 2
  expected <- meat * n / (n - p)
  result <- finite_sample_correction(meat, n, p, adjustment = "HC1")
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("finite_sample_correction() preserves dimensions", {
  meat <- diag(3)
  result <- finite_sample_correction(meat, n = 20, p = 3, adjustment = "HC1")
  expect_equal(dim(result), c(3L, 3L))
})

test_that("finite_sample_correction() errors when n <= p with HC1", {
  meat <- diag(3)
  expect_error(finite_sample_correction(meat, n = 3, p = 3, adjustment = "HC1"))
  expect_error(finite_sample_correction(meat, n = 2, p = 3, adjustment = "HC1"))
})

test_that("finite_sample_correction() HC1 approaches 1 for large n", {
  meat <- matrix(c(4, 1, 1, 9), nrow = 2, ncol = 2)
  result <- finite_sample_correction(meat, n = 10000, p = 2, adjustment = "HC1")
  expect_equal(result, meat, tolerance = 1e-3)
})

test_that("finite_sample_correction() works for 1x1 meat", {
  meat <- matrix(5, nrow = 1, ncol = 1)
  expected <- matrix(5 * 10 / 9, nrow = 1, ncol = 1)
  result <- finite_sample_correction(meat, n = 10, p = 1, adjustment = "HC1")
  expect_equal(result, expected, tolerance = 1e-6)
})

# ---- compute_sandwich() integration -----------------------------------------

test_that("compute_sandwich() matches Python for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  # Estimating equation for the mean: psi_i(theta) = y_i - theta
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  result <- compute_sandwich(psi, theta = ref$theta)

  # Python asymptotic_variance
  expect_equal(result[1, 1], ref$asymptotic_variance[[1]], tolerance = 1e-4)
})

test_that("compute_sandwich() matches Python for mean+variance EE", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  result <- compute_sandwich(psi, theta = ref$theta)

  # Check diagonal matches Python asymptotic variance
  expect_equal(dim(result), c(2L, 2L))
  expect_equal(diag(result), diag(ref$asymptotic_variance), tolerance = 1e-4)
})

# ---- compute_sandwich() as a top-level export -------------------------------
#
# compute_sandwich is one of the top-level entry points of the package, callable
# directly by a user who has already solved for theta and wants the sandwich
# without running MEstimator. The assertions below are of two kinds. The
# membership check on the deli namespace exports (bd-gxa4.5) gates the export
# itself: it holds only once compute_sandwich is exported. The remaining tests
# pin the user-facing contract the export must satisfy. compute_sandwich returns
# the asymptotic sandwich variance: the bread and meat are each scaled by n, so
# the returned matrix is not divided by n. That matrix equals the
# asymptotic_variance an MEstimator reports for the same estimating equation and
# root, dividing it by n recovers the standard-error-scale variance, and the
# numbers match Python delicatessen across every derivative method the two
# implementations share.

test_that("compute_sandwich is exported from the deli namespace", {
  expect_true("compute_sandwich" %in% getNamespaceExports("deli"))
})

test_that("compute_sandwich returns the MEstimator asymptotic variance", {
  ref <- load_fixture("sandwich_mean_variance")
  y <- ref$y

  psi <- function(theta) ee_mean_variance(theta, y = y)

  m <- estimate(MEstimator(psi, init = ref$init))
  sandwich <- compute_sandwich(psi, theta = m@theta, deriv_method = "capprox")

  # compute_sandwich returns the asymptotic (un-scaled by n) sandwich, which is
  # exactly what the estimator stores as its asymptotic variance.
  expect_equal(dim(sandwich), c(2L, 2L))
  expect_equal(
    unname(sandwich),
    unname(m@asymptotic_variance),
    tolerance = 1e-8
  )

  # Scaling the returned matrix by 1/n gives the standard-error-scale variance.
  expect_equal(unname(sandwich / m@n_obs), unname(m@variance), tolerance = 1e-8)
})

test_that("compute_sandwich matches Python for the mean-variance equation", {
  ref <- load_fixture("sandwich_mean_variance")
  y <- ref$y

  psi <- function(theta) ee_mean_variance(theta, y = y)

  for (method in c("exact", "capprox", "fapprox", "bapprox")) {
    sandwich <- compute_sandwich(psi, theta = ref$theta, deriv_method = method)
    expect_equal(
      unname(sandwich),
      as.matrix(ref$sandwich[[method]]),
      tolerance = 1e-6,
      label = paste0("sandwich_mean_variance: ", method)
    )
  }
})

test_that("compute_sandwich matches Python for a linear regression equation", {
  ref <- load_fixture("sandwich_regression")
  X <- as.matrix(ref$X)
  y <- ref$y

  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")

  # The finite-difference methods are compared here. The exact-autodiff path for
  # ee_regression is exercised by the mean-variance test above and tracked
  # separately in bd-3h2h; ee_regression does not yet preserve tangents through
  # its matrix operations, so it is left out of this loop.
  for (method in c("capprox", "fapprox", "bapprox")) {
    sandwich <- compute_sandwich(psi, theta = ref$theta, deriv_method = method)
    expect_equal(dim(sandwich), c(3L, 3L))
    expect_equal(
      unname(sandwich),
      as.matrix(ref$sandwich[[method]]),
      tolerance = 1e-6,
      label = paste0("sandwich_regression: ", method)
    )
  }
})

# deriv_method is case-insensitive --------------------------------------------
#
# Python lowercases every deriv_method comparison and accepts any case. deli
# compared case-sensitively, so "Exact" or "CAPPROX" raised "not supported".
# These tests pin case-insensitive acceptance and a clear error, listing all
# supported options, for a genuinely unknown method.

test_that("compute_sandwich() accepts deriv_method in any case", {
  ref <- load_fixture("sandwich_mean_variance")
  y <- ref$y
  psi <- function(theta) ee_mean_variance(theta, y = y)

  lower <- compute_sandwich(psi, theta = ref$theta, deriv_method = "capprox")
  upper <- compute_sandwich(psi, theta = ref$theta, deriv_method = "CAPPROX")
  mixed <- compute_sandwich(psi, theta = ref$theta, deriv_method = "CApprox")
  expect_equal(upper, lower)
  expect_equal(mixed, lower)
})

test_that("compute_sandwich() accepts an upper-case exact deriv_method", {
  ref <- load_fixture("sandwich_mean_variance")
  y <- ref$y
  psi <- function(theta) ee_mean_variance(theta, y = y)

  lower <- compute_sandwich(psi, theta = ref$theta, deriv_method = "exact")
  upper <- compute_sandwich(psi, theta = ref$theta, deriv_method = "Exact")
  expect_equal(upper, lower)
})

test_that("compute_sandwich() error for an unknown deriv_method lists exact", {
  ref <- load_fixture("sandwich_mean_variance")
  y <- ref$y
  psi <- function(theta) ee_mean_variance(theta, y = y)

  expect_error(
    compute_sandwich(psi, theta = ref$theta, deriv_method = "nonsense"),
    regexp = "exact"
  )
})

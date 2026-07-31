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

# Every condition this package raises goes through cli, so this one is an rlang
# warning rather than a base one, and cli::cli_warn() records no call because no
# site here passes one through to rlang::warn(). R/conditions.R keys the
# de-duplication scope on the class and the message, and relies on the absence
# of the call, so both are pinned here.
test_that("compute_bread() warns through cli when the bread contains NA", {
  psi <- function(theta) matrix(rep(NA_real_, 3) * theta[1], nrow = 1)

  w <- expect_warning(compute_bread(psi, theta = 1), "bread matrix contains NA")

  expect_s3_class(w, "rlang_warning")
  expect_null(conditionCall(w))
  expect_match(conditionMessage(w), "variance will not be calculated")
})

test_that("the NA bread warning carries the class the catalog documents", {
  # `deli_bread_na` is documented on `?deli-conditions` as the class a caller
  # matches to answer an NA bread, and every assertion about it read the prose
  # instead, which is the part a message is free to change. The class and the
  # wording are pinned together here.
  psi <- function(theta) matrix(rep(NA_real_, 3) * theta[1], nrow = 1)

  w <- expect_warning(
    compute_bread(psi, theta = 1),
    class = "deli_bread_na"
  )
  expect_equal(class(w)[[1L]], "deli_bread_na")
  flat <- gsub("\\s+", " ", conditionMessage(w))
  expect_match(flat, "bread matrix contains NA values", fixed = TRUE)
  expect_match(flat, "cannot be inverted", fixed = TRUE)
  expect_match(flat, "variance will not be calculated", fixed = TRUE)

  # `compute_sandwich()` converts an NA bread into the error the catalog names,
  # so a caller who cannot carry on without a variance matches that one. The
  # return has to be finite where it is asked for and NA only at a perturbed
  # point, since a return that is non-finite at `theta` is refused ahead of the
  # bread.
  at_theta <- function(theta) {
    matrix(rep(if (theta[1] == 0) 0 else NA_real_, 3), nrow = 1)
  }
  expect_error(
    compute_sandwich(at_theta, theta = 0),
    class = "deli_bread_not_invertible"
  )
})

# The finite-difference step is an absolute perturbation, so a large enough
# parameter magnitude drives it below the spacing of the doubles around `theta`
# and the Jacobian degrades, then vanishes. Nothing about the returned bread
# says so: `anyNA()` does not fire, `build_sandwich()` pseudo-inverts the zero
# matrix, and the reported standard error is exactly zero.
test_that("the bread survives a parameter magnitude that costs the step", {
  for (scale in c(5e6, 5e8, 5e9, 5e10)) {
    set.seed(42)
    y <- stats::rnorm(200, mean = scale, sd = scale / 5)
    psi <- function(theta) matrix(y - theta[1], nrow = 1)
    expected_variance <- sum((y - mean(y))^2) / 200^2

    for (method in c("capprox", "fapprox", "bapprox", "exact")) {
      expect_no_warning({
        m <- estimate(
          MEstimator(stacked_equations = psi, init = scale),
          deriv_method = method
        )
      })
      lab <- paste0(method, " at scale ", scale)
      expect_equal(unname(m@bread[1, 1]), 1, tolerance = 1e-3, label = lab)
      expect_equal(
        unname(m@variance[1, 1]),
        expected_variance,
        tolerance = 1e-3,
        label = lab
      )
    }
  }
})

test_that("a fit whose parameters need different steps matches exact", {
  # A design whose columns sit on different scales puts the three coefficients
  # at magnitudes 2, 2e6, and 5e3, so no two of them take the same step: the
  # first keeps `dx` itself, and the floor raises the other two to 4.3e-6 and
  # 1.1e-8. Each column of the bread therefore has to be divided by its own
  # parameter's step, and a bread built with any other correspondence is wrong
  # by up to the ratio between them.
  set.seed(11)
  n <- 250
  z <- stats::runif(n, 1, 3) / 1e6
  v <- stats::runif(n, 1, 3) / 1e3
  X <- cbind(1, z, v)
  y <- 2 + 2e6 * z + 5e3 * v + stats::rnorm(n)
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")

  fit <- function(method) {
    estimate(
      MEstimator(stacked_equations = psi, init = c(2, 2e6, 5e3)),
      deriv_method = method
    )
  }
  exact <- fit("exact")

  # The premise: the fit really does straddle the magnitude where the floor
  # engages, so the steps are not all `dx` and not all the floor.
  floor_threshold <- 1e-9 / (1e4 * .Machine$double.eps)
  expect_lt(abs(coef(exact)[[1]]), floor_threshold)
  expect_gt(min(abs(coef(exact)[-1])), floor_threshold)

  for (method in c("capprox", "fapprox", "bapprox")) {
    approx <- fit(method)
    # Only the bread depends on `deriv_method`, so the estimates must be bit
    # identical and every difference below belongs to the derivative.
    expect_identical(coef(approx), coef(exact))
    # Relative, entry by entry: the bread spans twelve orders of magnitude
    # here, so an absolute comparison would leave its small entries free. The
    # measured worst cases are 3.6e-05 for the bread and 5.6e-05 for the
    # standard errors.
    expect_lt(
      max(abs(unname(approx@bread / exact@bread) - 1)),
      1e-3,
      label = paste0("worst relative bread error, ", method)
    )
    expect_lt(
      max(abs(sqrt(diag(approx@variance)) / sqrt(diag(exact@variance)) - 1)),
      1e-3,
      label = paste0("worst relative standard error, ", method)
    )
  }
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
# membership check on the deli namespace exports gates the export
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
  # separately; ee_regression does not yet preserve tangents through
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

# ---- the per-equation list classification under the exact pass ---------------
#
# An estimating function may return one value per equation in a plain list, which
# an `lapply()` or `sapply()` over equations produces, and `summed_ee()` reduces
# that with `lapply(ef, sum)`: every tangent surface has a `sum()` method, so the
# reduction is a sum across observations within each equation whatever class the
# elements hold. Reading the class of the first element alone decided that on
# behalf of the whole list, so a list whose tangent-carrying element sat anywhere
# else fell through to the guard below and was refused a reduction the branch
# above would have carried out. The guard reads every element, and the
# classification reads every element too.
#
# The guard's own refusals stand. A list holding no scalar pair at all is the
# shape the abort describes, and a list holding no tangents anywhere is a lost
# derivative rather than an unsupported shape.

test_that("compute_bread() sums a per-equation list that leads with a constant", {
  # The first equation does not involve theta, so its element arrives as a plain
  # numeric vector with no tangent to read and a derivative of zero.
  y <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(rep(0, length(y)), theta[2] - y)
  # The same equations returned as one tangent-carrying matrix, differentiated by
  # finite differences, which is the reference the list shape has to reproduce.
  stacked <- function(theta) rbind(rep(0, length(y)), theta[2] - y)

  expect_equal(
    compute_bread(per_equation, c(2.5, 3.75), deriv_method = "exact"),
    compute_bread(stacked, c(2.5, 3.75), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

test_that("compute_bread() sums a per-equation list that leads with an array", {
  # c() on a scalar pair returns a tangent array, so an equation written that way
  # puts a PrimalTangentArray first and the scalar pair second.
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(c(theta[1] - y1), theta[2] - y2)
  stacked <- function(theta) rbind(theta[1] - y1, theta[2] - y2)

  expect_equal(
    compute_bread(per_equation, c(2.5, 3.75), deriv_method = "exact"),
    compute_bread(stacked, c(2.5, 3.75), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

test_that("compute_bread() still sums a per-equation list of scalar pairs", {
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(theta[1] - y1, theta[2] - y2)
  stacked <- function(theta) rbind(theta[1] - y1, theta[2] - y2)

  expect_equal(
    compute_bread(per_equation, c(2.5, 3.75), deriv_method = "exact"),
    compute_bread(stacked, c(2.5, 3.75), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

test_that("compute_bread() keeps the guard for the lists it cannot reduce", {
  # The boundary of the classification. A list holding no scalar pair keeps its
  # derivatives and reports the container shape; a list holding no tangents at
  # all reports the lost derivative. Neither is drawn into the reduction by
  # reading every element rather than the first.
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  arrays <- function(theta) list(c(theta[1] - y1), c(theta[2] - y2))
  expect_error(
    compute_bread(arrays, c(2.5, 3.75), deriv_method = "exact"),
    class = "deli_exact_unsupported_shape"
  )

  stripped <- function(theta) base::rbind(y1 - theta[1])
  expect_error(
    compute_bread(stripped, 2.5, deriv_method = "exact"),
    class = "deli_exact_tangent_lost"
  )
})

# ---- one element, one equation ----------------------------------------------
#
# Each element of a per-equation list reduces with a single `sum()`, so each
# contributes exactly one row to the bread. An element holding a block of
# several equations still reduces to one value, which means it silently stands
# in for all of the equations it holds and the bread comes back with fewer rows
# than the system has. Counting the elements against the parameters is what
# separates that from a list a reduction is correct for.

test_that("a list element holding a block of equations is refused", {
  # The mixed shape: a 2-by-n matrix block for the first two parameters beside a
  # scalar-pair equation for the third. Reducing it produced a 2-by-3 bread for
  # a 3-parameter system, which nothing downstream can tell from a correct one.
  y <- c(1, 2, 3)
  mixed <- function(theta) {
    list(rbind(theta[1] - y, theta[2] - y), theta[3] - y)
  }
  expect_error(
    compute_bread(mixed, c(1, 2, 3), deriv_method = "exact"),
    class = "deli_exact_unsupported_shape"
  )
  expect_error(
    compute_bread(mixed, c(1, 2, 3), deriv_method = "capprox"),
    class = "deli_exact_unsupported_shape"
  )
})

test_that("the element-count abort says what a list element stands for", {
  y <- c(1, 2, 3)
  mixed <- function(theta) {
    list(rbind(theta[1] - y, theta[2] - y), theta[3] - y)
  }
  err <- expect_error(
    compute_bread(mixed, c(1, 2, 3), deriv_method = "exact"),
    class = "deli_exact_unsupported_shape"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "one estimating equation", fixed = TRUE)
  expect_match(flat, "2 elements", fixed = TRUE)
})

# ---- one sample, one length -------------------------------------------------
#
# The equations of one system are evaluated at one sample, so a per-equation
# list holds one length. An element of another length is an equation built from
# the wrong observations, and the reduction hides it: each element is summed
# whatever its length, so the bread comes back with the documented shape and one
# row computed from a different number of contributions than the rest.

test_that("compute_bread() refuses a per-equation list of unequal lengths", {
  ragged <- function(theta) list(theta[1] - c(1, 2, 3, 4), theta[2] - c(2, 3))

  for (method in c("exact", "capprox", "fapprox", "bapprox")) {
    expect_error(
      compute_bread(ragged, c(2.5, 2.5), deriv_method = method),
      class = "deli_exact_unsupported_shape"
    )
  }
})

test_that("the length rule catches a block whose element count happens to fit", {
  # Two elements for two parameters passes the count, so the count alone let a
  # 2-by-n block stand beside a single equation and reduced the pair to a
  # 2-row bread for a system holding three equations. The block is twice the
  # width of its neighbor, which is what the length rule reads.
  y <- c(1, 2, 3)
  blocked <- function(theta) {
    list(rbind(theta[1] - y, theta[2] - y), theta[1] - y)
  }

  expect_error(
    compute_bread(blocked, c(1, 2), deriv_method = "capprox"),
    class = "deli_exact_unsupported_shape"
  )
})

test_that("the unequal-length abort says each equation holds one value per observation", {
  ragged <- function(theta) list(theta[1] - c(1, 2, 3, 4), theta[2] - c(2, 3))
  err <- expect_error(
    compute_bread(ragged, c(2.5, 2.5), deriv_method = "capprox"),
    class = "deli_exact_unsupported_shape"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "one value per observation", fixed = TRUE)
  expect_match(flat, "2 different lengths", fixed = TRUE)
})

test_that("an equal-length per-equation list is untouched by the length rule", {
  # The rule must not reach the shapes the reduction is correct for, including
  # the theta-free equation that arrives as a plain numeric vector of the same
  # length as the rest.
  y <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(rep(0, length(y)), theta[2] - y)
  stacked <- function(theta) rbind(rep(0, length(y)), theta[2] - y)

  expect_equal(
    compute_bread(per_equation, c(2.5, 3.75), deriv_method = "exact"),
    compute_bread(stacked, c(2.5, 3.75), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

# ---- the per-equation list under the finite-difference methods ---------------
#
# The package's rule is that nothing succeeds under the exact pass that would
# fail without it. A per-equation list evaluated at plain numbers is the same
# return the exact pass reduces, so the finite-difference methods reduce it the
# same way rather than handing `sum()` a list and reporting base R's
# `invalid 'type' (list) of argument`.

test_that("a per-equation list reduces under every finite-difference method", {
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(theta[1] - y1, theta[2] - y2)
  stacked <- function(theta) rbind(theta[1] - y1, theta[2] - y2)

  exact <- compute_bread(per_equation, c(2.5, 3.75), deriv_method = "exact")
  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_equal(
      compute_bread(per_equation, c(2.5, 3.75), deriv_method = method),
      exact,
      tolerance = 1e-6
    )
    expect_equal(
      compute_bread(per_equation, c(2.5, 3.75), deriv_method = method),
      compute_bread(stacked, c(2.5, 3.75), deriv_method = method),
      tolerance = 1e-6
    )
  }
})

# ---- the bread carries no names, whichever seam it came through --------------
#
# `estimate()` reads parameter names off the plain numeric evaluation it makes
# at the solved values, never off a differentiated one, so every intermediate
# the bread is assembled from is stripped rather than relying on the single
# `as.vector()` that reads the Jacobian column.

test_that("a psi that names its rows still produces an unnamed bread", {
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  named_rows <- function(theta) {
    out <- rbind(theta[1] - y1, theta[2] - y2)
    rownames(out) <- c("first", "second")
    out
  }

  for (method in c("exact", "capprox")) {
    bread <- compute_bread(named_rows, c(2.5, 3.75), deriv_method = method)
    expect_null(dimnames(bread))
  }
})

test_that("a psi returning a named per-equation list produces an unnamed bread", {
  # The list branch carried the element names into the reduction, and the
  # Jacobian column read off it kept them, so the bread came back with row
  # labels the numeric pass never reports.
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  named_list <- function(theta) {
    list(first = theta[1] - y1, second = theta[2] - y2)
  }

  for (method in c("exact", "capprox")) {
    bread <- compute_bread(named_list, c(2.5, 3.75), deriv_method = method)
    expect_null(dimnames(bread))
  }
})

# ---- the dimensionless tangent array under the exact pass -------------------
#
# `c()` on tangent-carrying values returns a tangent array whose primal has no
# `dim`, and summing that shape adds every element into one value rather than
# within each row. Adding across the whole return is the correct reduction for
# it: a dimensionless return is one estimating equation observed n times, which
# is the only reading `estimate()` gives it. The shape needs assembling with
# `c()` to arrive at all. A one-parameter estimating function that subtracts
# theta from a data vector returns a scalar pair instead, whose primal carries
# the vector and which is summed a step earlier.

test_that("an exact fit whose psi is assembled with c() matches its capprox twin", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    c(theta[1] - y[1], theta[1] - y[2], theta[1] - y[3])
  }

  exact <- m_estimate(stacked_equations = psi, init = 0, deriv_method = "exact")
  capprox <- m_estimate(
    stacked_equations = psi,
    init = 0,
    deriv_method = "capprox"
  )

  # Three contributions of one parameter, whose root is their mean.
  expect_equal(unname(exact@theta), mean(y[1:3]), tolerance = 1e-8)
  expect_equal(exact@theta, capprox@theta, tolerance = 1e-8)
  expect_equal(exact@bread, capprox@bread, tolerance = 1e-6)
  expect_equal(exact@variance, capprox@variance, tolerance = 1e-6)
})

test_that("compute_bread() sums a dimensionless tangent array across observations", {
  y <- c(1, 2, 4)
  psi <- function(theta) c(theta[1] - y[1], theta[1] - y[2], theta[1] - y[3])

  # Each contribution has derivative 1 with respect to theta, so the summed
  # Jacobian is 3 and the bread is its negation.
  expect_equal(
    compute_bread(psi, theta = mean(y), deriv_method = "exact"),
    matrix(-3)
  )
  expect_equal(
    compute_bread(psi, theta = mean(y), deriv_method = "exact"),
    compute_bread(psi, theta = mean(y), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

test_that("compute_sandwich() differentiates a c()-assembled psi exactly", {
  y <- c(1, 2, 4)
  psi <- function(theta) c(theta[1] - y[1], theta[1] - y[2], theta[1] - y[3])

  expect_equal(
    compute_sandwich(psi, theta = mean(y), deriv_method = "exact"),
    compute_sandwich(psi, theta = mean(y), deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

# ---- compute_sandwich() judges the estimating-function return ----------------
#
# compute_sandwich() assembles the sandwich at a point the caller states is the
# root, so it never runs the validation estimate() runs before solving. The
# returns that validation rejects reach the bread and the meat here instead, and
# each of them either fails deep in base R with a message about an argument the
# caller never named or, worse, returns a matrix of the documented shape built
# from a Jacobian that is not one. The same judgment belongs on this entry
# point, with one allowance: more estimating equations than parameters is the
# over-identified GMM system, whose variance a GMMEstimator reports by
# assembling exactly this sandwich from exactly this rectangular bread.

test_that("compute_sandwich() rejects fewer estimating equations than parameters", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  # One equation for two parameters. The bread is 1-by-2, and the pseudo-inverse
  # in build_sandwich() turns that into a 2-by-2 matrix carrying the documented
  # shape and none of the documented meaning.
  psi <- function(theta) y - theta[1] - theta[2]

  # The wording is left to the abort that judges it; what the count of equations
  # and the name of the argument that returned them have to survive is the move
  # from estimate() to this entry point.
  expect_error(
    compute_sandwich(psi, theta = c(0.5, 0.5)),
    regexp = "estimating equation",
    class = "deli_psi_shape_error"
  )
  expect_error(
    compute_sandwich(psi, theta = c(0.5, 0.5)),
    regexp = "stacked_equations"
  )
})

test_that("compute_sandwich() rejects a NULL return", {
  # Without the judgment this dies inside base::matrix(), reporting `data` as
  # the offending argument. The caller passed no `data`.
  psi <- function(theta) NULL

  expect_error(compute_sandwich(psi, theta = 0), regexp = "stacked_equations")
  expect_error(compute_sandwich(psi, theta = 0), regexp = "NULL")
})

test_that("compute_sandwich() rejects a non-numeric return", {
  # Without the judgment this dies inside rowSums(), reporting `x`.
  psi <- function(theta) matrix("a", nrow = 1, ncol = 3)

  expect_error(compute_sandwich(psi, theta = 0), regexp = "stacked_equations")
  expect_error(compute_sandwich(psi, theta = 0), regexp = "numeric")
})

test_that("compute_sandwich() rejects a non-finite return", {
  # An estimating function that is not finite at the point it is evaluated at
  # gives an all-NA bread and an all-NA meat, and the sandwich built from them
  # is reported as a bread that cannot be inverted. The bread is a symptom; the
  # return is the cause, and it is the one the caller can act on.
  psi <- function(theta) matrix(rep(NA_real_, 3) * theta[1], nrow = 1)

  expect_error(compute_sandwich(psi, theta = 1), regexp = "stacked_equations")
  expect_error(compute_sandwich(psi, theta = 1), regexp = "non-finite")
})

test_that("compute_sandwich() accepts an over-identified system", {
  skip_if_not_installed("MASS")
  # Three moment conditions for two parameters. The bread is 3-by-2, which
  # build_sandwich() pseudo-inverts into the identity-weighted GMM variance.
  # That is the same variance a GMMEstimator reports, from the same two
  # matrices, so rejecting the shape here would withdraw a result the package
  # already returns elsewhere.
  set.seed(20)
  n <- 150
  x <- stats::rnorm(n)
  y <- 1 + 0.6 * x + stats::rnorm(n)
  X <- cbind(1, x)

  psi <- function(theta) {
    r <- as.vector(y - X %*% theta)
    rbind(r, r * x, r * x^2)
  }

  fit <- gmm_estimate(stacked_equations = psi, init = c(0, 0))
  sandwich <- compute_sandwich(psi, theta = fit@theta)

  expect_equal(dim(sandwich), c(2L, 2L))
  expect_equal(
    sandwich,
    unname(fit@asymptotic_variance),
    tolerance = 1e-8
  )
})

test_that("compute_sandwich() accepts a one-parameter psi returning a bare vector", {
  # A dimensionless return is one estimating equation observed n times, which is
  # a legitimate shape for a one-parameter fit and must survive the judgment
  # that rejects it for two parameters.
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  bare <- compute_sandwich(function(theta) y - theta[1], theta = mean(y))
  shaped <- compute_sandwich(
    function(theta) matrix(y - theta[1], nrow = 1),
    theta = mean(y)
  )

  expect_equal(dim(bare), c(1L, 1L))
  expect_equal(bare, shaped, tolerance = 1e-8)
  expect_equal(
    bare[1, 1],
    sum((y - mean(y))^2) / length(y),
    tolerance = 1e-6
  )
})

# ---- the per-equation list through compute_sandwich() ------------------------
#
# A list of one element per equation is the return compute_bread() reduces, and
# the two entry points that read a variance at a point supplied as the root are
# meant to take the same returns. compute_sandwich() judged the raw list and
# refused it as a non-numeric return before the bread was ever built, so a
# caller who wrote an estimating function for the documented shape could reach
# the Jacobian and not the sandwich.
#
# The bread reduces each element with one `sum()`, which is the derivative it
# needs; the meat needs the same equations unreduced, so the elements are bound
# into the p-by-n matrix here. Both halves accept exactly the lists the bread's
# element-count rule accepts.

test_that("compute_sandwich() assembles a sandwich from a per-equation list", {
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(theta[1] - y1, theta[2] - y2)
  stacked <- function(theta) rbind(theta[1] - y1, theta[2] - y2)
  theta <- c(mean(y1), mean(y2))

  expect_equal(
    compute_sandwich(per_equation, theta = theta),
    compute_sandwich(stacked, theta = theta),
    tolerance = 1e-8
  )
})

test_that("a per-equation list reaches the sandwich under every deriv_method", {
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3, 4, 6)
  per_equation <- function(theta) list(theta[1] - y1, theta[2] - y2)
  stacked <- function(theta) rbind(theta[1] - y1, theta[2] - y2)
  theta <- c(mean(y1), mean(y2))

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_equal(
      compute_sandwich(per_equation, theta = theta, deriv_method = method),
      compute_sandwich(stacked, theta = theta, deriv_method = method),
      tolerance = 1e-8
    )
  }
  # The stacked reference is differentiated by finite differences, because its
  # `rbind()` is the base one when the estimating function is written outside
  # the package and so carries no tangents through the exact pass.
  expect_equal(
    compute_sandwich(per_equation, theta = theta, deriv_method = "exact"),
    compute_sandwich(stacked, theta = theta, deriv_method = "capprox"),
    tolerance = 1e-6
  )
})

test_that("compute_sandwich() refuses the list shapes the bread refuses", {
  # An element holding a block of several equations sums to a single value, so
  # it stands in for every equation it holds. The bread rejects it by counting
  # elements against parameters, and the sandwich rejects it on the same count
  # rather than building a meat the bread cannot be paired with.
  y <- c(1, 2, 3)
  mixed <- function(theta) {
    list(rbind(theta[1] - y, theta[2] - y), theta[3] - y)
  }

  expect_error(
    compute_sandwich(mixed, theta = c(1, 2, 3)),
    class = "deli_exact_unsupported_shape"
  )
  expect_error(
    compute_sandwich(mixed, theta = c(1, 2, 3), deriv_method = "exact"),
    class = "deli_exact_unsupported_shape"
  )
})

test_that("compute_sandwich() refuses a per-equation list of unequal lengths", {
  # `rbind()` recycles a short element up to the width of the longest one, and a
  # length that divides that width is recycled with nothing reported, so the
  # meat would be built from repeated contributions and carry the shape of a
  # covariance matrix without the meaning. The bread refuses the same list on
  # the same rule, so both halves report one class.
  y1 <- c(1, 2, 3, 4)
  y2 <- c(2, 3)
  ragged <- function(theta) list(theta[1] - y1, theta[2] - y2)

  err <- expect_error(
    compute_sandwich(ragged, theta = c(2.5, 2.5)),
    class = "deli_exact_unsupported_shape"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "one value per observation", fixed = TRUE)
})

# ---- the frame those aborts report -------------------------------------------
#
# The judgment above is made in a helper, so each of its aborts reported a frame
# naming that helper's own parameters rather than the call the caller made. The
# entry point the caller can act on is compute_sandwich(), and it is the one the
# report names.

reported_call <- function(err) {
  paste(deparse(conditionCall(err)), collapse = " ")
}

test_that("the NULL-return abort names compute_sandwich()", {
  psi <- function(theta) NULL
  # These aborts share a call and a class, so the wording is what tells them
  # apart and the frame assertion needs it to know which one it caught.
  err <- expect_error(compute_sandwich(psi, theta = 0), regexp = "NULL")

  expect_match(reported_call(err), "compute_sandwich(", fixed = TRUE)
  expect_false(grepl("check_psi_at_theta", reported_call(err), fixed = TRUE))
})

test_that("the shortfall abort names compute_sandwich()", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) y - theta[1] - theta[2]
  err <- expect_error(
    compute_sandwich(psi, theta = c(0.5, 0.5)),
    class = "deli_psi_shape_error"
  )

  expect_match(reported_call(err), "compute_sandwich(", fixed = TRUE)
  expect_false(grepl("check_psi_at_theta", reported_call(err), fixed = TRUE))
})

test_that("the non-finite abort names compute_sandwich()", {
  psi <- function(theta) matrix(rep(NA_real_, 3) * theta[1], nrow = 1)
  err <- expect_error(compute_sandwich(psi, theta = 1), regexp = "non-finite")

  expect_match(reported_call(err), "compute_sandwich(", fixed = TRUE)
  expect_false(grepl("check_psi_at_theta", reported_call(err), fixed = TRUE))
})

# ---- the bread that cannot be inverted ---------------------------------------
#
# compute_sandwich() returns a covariance matrix or it fails. Three breads
# cannot be turned into one, and each of them used to be reported as something
# other than a failure of this call: an NA bread warned and returned `NULL`, so
# the caller indexed into `NULL` and got `NULL` back one subscript at a time; a
# rank-deficient bread under `allow_pinv = FALSE` cleared base R's inversion
# tolerance and returned finite standard errors that mean nothing; and a
# rectangular bread under the same setting died inside base::solve() against an
# argument named `a`. All three now carry `deli_bread_not_invertible`.
#
# estimate() is the other caller of the same assembly and keeps its own
# contract: a fit whose bread holds NA warns and comes back with no variance,
# which is a state the accessors name. Only the entry point that has nothing but
# the matrix to return fails.

test_that("compute_sandwich() rejects an NA bread with one classed error", {
  # Finite where the meat reads it and not finite where the difference quotient
  # does, so the bread is the only part of the sandwich holding an NA.
  psi <- function(theta) {
    matrix(rep(if (theta[1] == 1) 1 else NA_real_, 3), nrow = 1)
  }

  err <- expect_error(
    compute_sandwich(psi, theta = 1),
    class = "deli_bread_not_invertible"
  )

  expect_match(conditionMessage(err), "bread")
  expect_match(reported_call(err), "compute_sandwich(", fixed = TRUE)
})

test_that("the NA bread is reported once rather than warned about as well", {
  psi <- function(theta) {
    matrix(rep(if (theta[1] == 1) 1 else NA_real_, 3), nrow = 1)
  }

  caught <- collect_warnings(
    expect_error(
      compute_sandwich(psi, theta = 1),
      class = "deli_bread_not_invertible"
    )
  )

  expect_length(caught, 0)
})

test_that("estimate() still reports an NA bread as a fit with no variance", {
  # Zero at the starting values, so the point is a root and the bread is the
  # only thing the fit has to report on.
  psi <- function(theta) {
    matrix(rep(if (theta[1] == 0) 0 else NA_real_, 3), nrow = 1)
  }

  expect_warning(
    fit <- estimate(
      MEstimator(stacked_equations = psi, init = 0),
      solver = function(stacked_equations, init) init
    ),
    "bread matrix contains NA"
  )

  expect_null(fit@variance)
})

# ---- the rank reading is the caller's, not the sandwich's ---------------------
#
# `allow_pinv = FALSE` says a bread that has no inverse is to be refused rather
# than pseudo-inverted, and what has no inverse is decided by the solve. A rank
# read ahead of it, at the 1e-7 relative tolerance `qr()` applies, refuses far
# more than that: a bread whose condition number runs past roughly 1e7 has an
# inverse, and base R hands it over, and a caller applying a more discriminating
# gate of their own never reached it.

# A design whose fourth column is the sum of two others is rank deficient
# analytically, and the bread of a logistic fit on it is read as rank deficient
# too. The finite differences do not reproduce the dependence exactly, though:
# the round-off in the difference quotient perturbs the dependent row enough
# that `solve()` inverts the matrix without complaint.
collinear_sandwich_case <- function() {
  set.seed(1)
  n <- 200
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, 1 / (1 + exp(-(0.5 + 0.8 * x1))))
  X <- cbind(1, x1, x2, x1 + x2)
  list(
    psi = function(theta) {
      ee_regression(theta, X = X, y = y, model = "logistic")
    },
    theta = c(0.5, 0.4, 0, 0)
  )
}

test_that("a bread solve() inverts is not refused for its rank reading", {
  case <- collinear_sandwich_case()
  bread <- compute_bread(case$psi, case$theta) / 200

  # The premise: `qr()` reads this bread as rank deficient at its own tolerance,
  # and base R inverts it without complaint.
  expect_lt(qr(bread)$rank, ncol(bread))
  expect_no_error(solve(bread))

  expect_no_error(
    compute_sandwich(case$psi, theta = case$theta, allow_pinv = FALSE)
  )
})

test_that("the same bread reaches no pseudo-inverse when one is allowed", {
  case <- collinear_sandwich_case()

  # `allow_pinv = TRUE` reaches `MASS::ginv()` only where `solve()` fails, and
  # it does not fail here, so the setting makes no difference to what comes
  # back.
  expect_identical(
    compute_sandwich(case$psi, theta = case$theta, allow_pinv = FALSE),
    compute_sandwich(case$psi, theta = case$theta)
  )
})

# A stacked system whose reported block is well conditioned and whose nuisance
# block carries the dependent column above. The two blocks share no parameters,
# so the bread is block diagonal and the reported block of the sandwich is the
# sandwich of the reported equations on their own. That is the shape the rank
# pre-check was measured costing: the whole matrix reads as rank deficient, the
# effects the caller reports are identified to every digit they had, and the
# refusal named a direction nobody was asking about.
nuisance_collinear_case <- function() {
  set.seed(4)
  n <- 300
  w <- stats::rnorm(n, 3, 2)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, 1 / (1 + exp(-(0.5 + 0.8 * x1))))
  X <- cbind(1, x1, x2, x1 + x2)
  list(
    n = n,
    psi = function(theta) {
      rbind(
        ee_mean_variance(theta[1:2], y = w),
        ee_regression(theta[3:6], X = X, y = y, model = "logistic")
      )
    },
    reported = function(theta) ee_mean_variance(theta, y = w),
    theta = c(mean(w), stats::var(w) * (n - 1) / n, 0.5, 0.4, 0, 0)
  )
}

test_that("a bread whose dependence sits in a nuisance block is inverted", {
  skip_if_not_installed("MASS")
  case <- nuisance_collinear_case()
  bread <- compute_bread(case$psi, case$theta) / case$n

  expect_lt(qr(bread)$rank, ncol(bread))
  expect_no_error(solve(bread))

  sandwich <- compute_sandwich(
    case$psi,
    theta = case$theta,
    allow_pinv = FALSE
  )
  reported <- sandwich[1:2, 1:2]

  # Against the pseudo-inverse the same call would have taken with
  # `allow_pinv = TRUE`, had the solve failed.
  meat <- compute_meat(case$psi(case$theta)) / case$n
  pinv <- MASS::ginv(bread)
  expect_equal(
    reported,
    (pinv %*% meat %*% t(pinv))[1:2, 1:2],
    tolerance = 1e-10
  )

  # And against the fit of the reported equations alone, which is the full-rank
  # system the nuisance block was stacked onto.
  expect_equal(
    reported,
    compute_sandwich(case$reported, theta = case$theta[1:2]),
    tolerance = 1e-10
  )
})

test_that("an exactly singular bread is still refused without the pinv", {
  # The second equation does not move when either parameter does, so its row of
  # the bread is zero however the difference quotient rounds, and `solve()`
  # fails on the matrix outright.
  psi <- function(theta) rbind(c(1, 2, 3) - theta[1], rep(0, 3))

  err <- expect_error(
    compute_sandwich(psi, theta = c(2, 0), allow_pinv = FALSE),
    class = "deli_bread_not_invertible"
  )

  expect_match(conditionMessage(err), "singular")
  expect_match(conditionMessage(err), "allow_pinv")
  # The rank reading is made on the failure path, where it says which
  # directions the matrix lost rather than whether to refuse at all.
  expect_match(conditionMessage(err), "rank is 1 of 2", fixed = TRUE)
})

test_that("a bread whose columns differ only in scale is refused by name", {
  # Columns that are independent to the factorization's tolerance and whose
  # scales differ by more than the reciprocal condition number `solve()`
  # accepts. Nothing is rank deficient here, so the refusal has no direction to
  # name and says what it does know.
  bread <- diag(c(1, 1e-300))
  meat <- diag(2)

  err <- expect_error(
    build_sandwich(bread, meat, allow_pinv = FALSE),
    class = "deli_bread_not_invertible"
  )

  expect_match(conditionMessage(err), "full rank", fixed = TRUE)
  expect_match(conditionMessage(err), "allow_pinv")
})

test_that("a rectangular bread is refused with the pseudo-inverse named", {
  set.seed(42)
  counts <- stats::rpois(200, 3)
  psi <- function(theta) {
    rbind(counts - theta[1], (counts - theta[1])^2 - theta[1])
  }

  err <- expect_error(
    compute_sandwich(psi, theta = 3, allow_pinv = FALSE),
    class = "deli_bread_not_invertible"
  )

  # base::solve() reported this as "'a' (2 x 1) must be square", naming an
  # argument of its own and saying nothing about the system that produced it.
  expect_match(conditionMessage(err), "over-identified")
  expect_match(conditionMessage(err), "allow_pinv")
  expect_false(grepl("must be square", conditionMessage(err), fixed = TRUE))
})

test_that("an over-identified fit still refuses the pseudo-inverse by class", {
  set.seed(42)
  counts <- stats::rpois(200, 3)
  psi <- function(theta) {
    rbind(counts - theta[1], (counts - theta[1])^2 - theta[1])
  }

  expect_error(
    gmm_estimate(stacked_equations = psi, init = 1, allow_pinv = FALSE),
    class = "deli_bread_not_invertible"
  )
})

# ---- the summed estimating equations the bread differentiates ----------------
#
# The bread is the Jacobian of the summed estimating equations, so every
# perturbed evaluation it makes is reduced to one value per equation the moment
# it is built. Deriving that reduction from the full p-by-n return costs the
# whole matrix 2p times over for arithmetic that is linear in it. A caller who
# can write the sums directly supplies them as `summed_equations`, and the meat
# still takes the one full evaluation it needs. The two must be the same system,
# which is what the agreement check below judges.

reducer_case <- function(n = 80) {
  set.seed(11)
  X <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  y <- as.vector(X %*% c(0.5, 1.5, -0.75)) + stats::rnorm(n)
  theta <- as.vector(qr.solve(X, y))
  list(
    X = X,
    y = y,
    theta = theta,
    psi = function(th) t(X * (y - as.vector(X %*% th))),
    # The row sums of that psi, written as one matrix product rather than as a
    # reduction over the p-by-n matrix. Under a finite-difference method the
    # argument is plain numeric, so the whole thing is `t(X) %*% resid`.
    summed = function(th) as.vector(t(X) %*% (y - as.vector(X %*% th)))
  )
}

test_that("a supplied reducer gives the sandwich the matrix route gives", {
  case <- reducer_case()

  expect_equal(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = case$summed
    ),
    compute_sandwich(case$psi, theta = case$theta),
    tolerance = 1e-5
  )
})

test_that("compute_bread() takes the same reducer", {
  case <- reducer_case()

  # The two reductions sum the same values in different orders, so their
  # difference quotients differ by rounding rather than exactly.
  expect_equal(
    compute_bread(case$psi, case$theta, summed_equations = case$summed),
    compute_bread(case$psi, case$theta),
    tolerance = 1e-5
  )
})

test_that("a supplied reducer is what the perturbed evaluations go through", {
  case <- reducer_case()
  calls <- 0L
  counted <- function(th) {
    calls <<- calls + 1L
    case$psi(th)
  }

  compute_sandwich(counted, theta = case$theta, summed_equations = case$summed)
  # The meat's one evaluation, and no other: the bread went through the
  # reduction it was handed.
  expect_identical(calls, 1L)

  calls <- 0L
  compute_sandwich(counted, theta = case$theta)
  # The meat's one evaluation, and two per parameter for the central
  # differences, which is what the argument exists to avoid.
  expect_identical(calls, 1L + 2L * length(case$theta))
})

test_that("the derived reduction is what an unsupplied reducer leaves in place", {
  case <- reducer_case()

  expect_identical(
    compute_sandwich(case$psi, theta = case$theta, summed_equations = NULL),
    compute_sandwich(case$psi, theta = case$theta)
  )
  expect_identical(
    compute_bread(case$psi, case$theta, summed_equations = NULL),
    compute_bread(case$psi, case$theta)
  )
})

test_that("a reducer that sums another system is refused", {
  case <- reducer_case()
  # The second equation carries a term the estimating function does not, so the
  # bread would be the Jacobian of one system and the meat the cross-product of
  # another, and the returned matrix would carry the shape of a covariance
  # without the meaning.
  wrong <- function(th) case$summed(th) + c(0, 5, 0)

  err <- expect_error(
    compute_sandwich(case$psi, theta = case$theta, summed_equations = wrong),
    class = "deli_summed_equations_disagree"
  )

  expect_match(conditionMessage(err), "summed_equations")
  expect_match(conditionMessage(err), "stacked_equations")
})

test_that("a reducer returning the wrong number of values is refused", {
  case <- reducer_case()
  short <- function(th) case$summed(th)[1:2]

  # A return of the wrong length is a shape the sandwich cannot be assembled
  # from rather than a disagreement between two systems, so it carries the
  # family class alone.
  err <- expect_error(
    compute_sandwich(case$psi, theta = case$theta, summed_equations = short),
    class = "deli_summed_equations_error"
  )

  expect_match(conditionMessage(err), "2")
  expect_match(conditionMessage(err), "3")
  expect_false(inherits(err, "deli_summed_equations_disagree"))
})

test_that("a reducer that is not a function is refused before it is called", {
  case <- reducer_case()

  err <- expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = c(1, 2, 3)
    ),
    class = "deli_summed_equations_error"
  )
  expect_match(conditionMessage(err), "summed_equations")

  expect_error(
    compute_bread(case$psi, case$theta, summed_equations = c(1, 2, 3)),
    class = "deli_summed_equations_error"
  )
})

test_that("a reducer returning something other than numbers is refused", {
  case <- reducer_case()
  as_text <- function(th) as.character(case$summed(th))
  as_list <- function(th) as.list(case$summed(th))

  for (bad in list(as_text, as_list)) {
    expect_error(
      compute_sandwich(case$psi, theta = case$theta, summed_equations = bad),
      class = "deli_summed_equations_error"
    )
    expect_error(
      compute_bread(case$psi, case$theta, summed_equations = bad),
      class = "deli_summed_equations_error"
    )
  }
})

test_that("compute_bread() refuses a reduction shorter than the parameters", {
  case <- reducer_case()
  short <- function(th) case$summed(th)[1:2]

  # It has no evaluation of the estimating functions to compare against, so the
  # shape is the whole of what it can read. Unread, this returned a two-by-three
  # bread, which has no inverse at any rank and says nothing about why.
  err <- expect_error(
    compute_bread(case$psi, case$theta, summed_equations = short),
    class = "deli_summed_equations_error"
  )
  expect_match(conditionMessage(err), "2")
  expect_match(conditionMessage(err), "3")
})

test_that("a reducer that returns NA is refused rather than differentiated", {
  case <- reducer_case()
  all_missing <- function(th) rep(NA_real_, 3)
  one_missing <- function(th) c(NA_real_, case$summed(th)[-1])

  # The whole return: this reported base R's `argument is of length zero`, from
  # the comparison itself, because there was no largest disagreement to name.
  err <- expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = all_missing
    ),
    class = "deli_summed_equations_disagree"
  )
  expect_match(conditionMessage(err), "NA")

  # One equation of it: this passed the comparison, which reads the largest
  # disagreement and ignores the missing one, and surfaced two steps later as a
  # bread that could not be inverted.
  partial <- expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = one_missing
    ),
    class = "deli_summed_equations_disagree"
  )
  expect_match(conditionMessage(partial), "equation 1", fixed = TRUE)
})

test_that("the narrower reducer class leads the family it belongs to", {
  case <- reducer_case()
  wrong <- function(th) case$summed(th) + c(0, 5, 0)

  err <- expect_error(
    compute_sandwich(case$psi, theta = case$theta, summed_equations = wrong),
    class = "deli_summed_equations_disagree"
  )
  expect_s3_class(err, "deli_summed_equations_error")
  expect_identical(class(err)[[1]], "deli_summed_equations_disagree")
})

test_that("the agreement check can be turned off", {
  case <- reducer_case()
  wrong <- function(th) case$summed(th) + c(0, 5, 0)

  expect_no_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = wrong,
      check_summed_equations = FALSE
    )
  )
})

test_that("a tangent-safe reducer differentiates exactly", {
  case <- reducer_case()
  # `t(X) %*% resid` carries tangents through the registered `%*%` and `t()`
  # methods, so the same reducer serves the exact pass. The crossprod shortcut
  # would not.
  expect_equal(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      deriv_method = "exact",
      summed_equations = case$summed
    ),
    compute_sandwich(case$psi, theta = case$theta, deriv_method = "exact")
  )
})

test_that("a reducer that rescales the same system is not what the check sees", {
  case <- reducer_case()
  # The comparison is made at `theta`, which the caller states is the root, so
  # the sums it compares are at rounding there and a reduction that is a
  # multiple of the right one agrees with it. The bread it yields is that
  # multiple of the right bread. The check sees a reduction of some other
  # system, not a rescaling of this one, and the documentation says so.
  rescaled <- function(th) case$summed(th) * 2

  expect_no_error(
    compute_sandwich(case$psi, theta = case$theta, summed_equations = rescaled)
  )
})

test_that("the argument and its shape are read with the comparison turned off", {
  case <- reducer_case()
  as_text <- function(th) as.character(case$summed(th))
  short <- function(th) case$summed(th)[1:2]

  # Neither reading is a comparison against the estimating functions, so both
  # stand where the comparison is off: the argument here, and the shape in
  # compute_bread(), which reads what it can know without evaluating them.
  expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = c(1, 2, 3),
      check_summed_equations = FALSE
    ),
    class = "deli_summed_equations_error"
  )
  expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = as_text,
      check_summed_equations = FALSE
    ),
    class = "deli_summed_equations_error"
  )
  expect_error(
    compute_sandwich(
      case$psi,
      theta = case$theta,
      summed_equations = short,
      check_summed_equations = FALSE
    ),
    class = "deli_summed_equations_error"
  )
})

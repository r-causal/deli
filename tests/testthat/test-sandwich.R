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
  err <- expect_error(compute_sandwich(psi, theta = 0))

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
  err <- expect_error(compute_sandwich(psi, theta = 1))

  expect_match(reported_call(err), "compute_sandwich(", fixed = TRUE)
  expect_false(grepl("check_psi_at_theta", reported_call(err), fixed = TRUE))
})

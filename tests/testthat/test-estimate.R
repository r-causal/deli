# ---- estimate() for mean EE --------------------------------------------------

test_that("estimate() finds correct theta for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
})

test_that("estimate() computes correct variance for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(unname(m@variance[1, 1]), ref$variance[[1]], tolerance = 1e-4)
})

test_that("estimate() computes correct asymptotic variance for mean EE", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(
    unname(m@asymptotic_variance[1, 1]),
    ref$asymptotic_variance[[1]],
    tolerance = 1e-4
  )
})

# ---- estimate() for mean+variance EE ----------------------------------------

test_that("estimate() finds correct theta for mean+variance EE", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
})

test_that("estimate() computes correct variance for mean+variance EE", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(unname(diag(m@variance)), diag(ref$variance), tolerance = 1e-4)
})

# ---- estimate() sets n_obs ---------------------------------------------------

test_that("estimate() sets n_obs correctly", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(m@n_obs, 5L)
})

# ---- estimate() sets bread and meat ------------------------------------------

test_that("estimate() populates bread and meat", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_true(!is.null(m@bread))
  expect_true(!is.null(m@meat))
  expect_true(is.matrix(m@bread))
  expect_true(is.matrix(m@meat))
})

# ---- estimate() returns MEstimator ------------------------------------------

test_that("estimate() returns an MEstimator object", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  result <- estimate(m)

  expect_s3_class(result, "deli::MEstimator")
})

# ---- estimate() validates psi(init) before solving --------------------------
#
# Malformed psi/init combinations previously failed deep inside the solver with
# opaque rootSolve internals, or silently returned the initial values. estimate()
# now evaluates the estimating functions once at the starting values and aborts
# with an informative message.

test_that("estimate() rejects a psi/init dimension mismatch", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0, 0))
  expect_error(estimate(m), "1 estimating equation.*2 parameter")
})

test_that("estimate() rejects a psi that is non-finite at init", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y * NaN - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m), "non-finite")
})

test_that("estimate() rejects a psi that returns NULL", {
  psi <- function(theta) NULL
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m), "NULL")
})

test_that("estimate() rejects non-finite data at init under the lm solver", {
  # Previously the lm solver silently returned theta = init with NULL variance.
  psi <- function(theta) matrix(c(1, 2, NA) - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m, solver = "lm"), "non-finite")
})

test_that("estimate() rejects a custom solver that returns the wrong shape", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  bad_solver <- function(stacked_equations, init) c(1, 2, 3)
  expect_error(estimate(m, solver = bad_solver), "numeric vector of length 1")
})

test_that("estimate() names the list return a fit has no support for", {
  # `compute_bread()` and `compute_sandwich()` both take a per-equation list, so
  # a caller who built one for them reasonably expects a fit to take it as well.
  # No fit does: the summed equations reach `sum()` on a dimensionless list,
  # which is base R's `invalid 'type' (list) of argument`, and the generic
  # refusal of a non-numeric return names the type without saying that the entry
  # points differ. Both estimator paths reach the same refusal.
  psi <- function(theta) {
    y <- c(1, 2, 3, 4, 5)
    list(y - theta[1], (y - theta[1])^2 - theta[2])
  }

  m <- expect_error(
    estimate(MEstimator(stacked_equations = psi, init = c(0, 1))),
    class = "deli_psi_list_unsupported"
  )
  g <- expect_error(
    estimate(GMMEstimator(stacked_equations = psi, init = c(0, 1))),
    class = "deli_psi_list_unsupported"
  )

  for (err in list(m, g)) {
    flat <- gsub("\\s+", " ", conditionMessage(err))
    expect_match(flat, "p-by-n", fixed = TRUE)
    expect_match(flat, "rbind", fixed = TRUE)
    expect_match(flat, "compute_bread", fixed = TRUE)
    expect_match(flat, "compute_sandwich", fixed = TRUE)
    # The family class leads with the narrower one, as every other refused
    # return does.
    expect_s3_class(err, "deli_psi_return_error")
  }
})

test_that("estimate() refuses a maxiter that is not a single positive count", {
  # The budget reached the solver unjudged, beside a `dx` and a `deriv_method`
  # that are both checked. A vector budget pluralized the non-convergence report
  # on its own length rather than on the budget, and a budget that is not a
  # number failed inside the solver against an argument the caller never wrote.
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = 0)
  g <- GMMEstimator(stacked_equations = psi, init = 0)

  budgets <- list(c(10, 20), "many", 0, -5, NA_integer_, 2.5, Inf, numeric(0))
  for (budget in budgets) {
    expect_error(estimate(m, maxiter = budget), "maxiter")
    expect_error(estimate(g, maxiter = budget), "maxiter")
  }
})

test_that("estimate() reports the budget it refuses against the caller", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = 0)

  err <- expect_error(estimate(m, maxiter = c(10, 20)))
  expect_match(
    reported_call(err),
    "estimate",
    fixed = TRUE
  )
  expect_no_error(estimate(m, maxiter = 100))
})

# ---- estimate(deriv_method = "exact") for built-in EEs ----------------------
#
# Exact-mode bread for a built-in estimating equation must agree with the
# central-difference (capprox) bread computed on the same fitted estimator. The
# cross-check uses tolerance 1e-5 because capprox carries roughly 1e-6 numerical
# error at the default step size (floating-point cancellation), consistent with
# the autodiff-vs-approx convention in test-autodiff.R; exact autodiff itself is
# near machine precision, so the tolerance is bounded by the approximation it is
# compared against. These EEs flow through the masked matrix()/rbind() surface
# without stripping tangents (no as.numeric() applied to a theta-carrying value),
# so they exercise the pt_flatten/bind normalization.

test_that("estimate() exact bread matches capprox for ee_mean", {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) ee_mean(theta, y = y)

  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  m_capprox <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "capprox"
  )

  expect_equal(m_exact@theta, m_capprox@theta, tolerance = 1e-9)
  expect_equal(m_exact@bread, m_capprox@bread, tolerance = 1e-5)
  expect_equal(m_exact@variance, m_capprox@variance, tolerance = 1e-5)
})

test_that("estimate() exact bread matches capprox for ee_mean_variance (rbind path)", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y
  psi <- function(theta) ee_mean_variance(theta, y = y)

  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  m_capprox <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "capprox"
  )

  expect_equal(m_exact@theta, m_capprox@theta, tolerance = 1e-9)
  expect_equal(m_exact@bread, m_capprox@bread, tolerance = 1e-5)
  expect_equal(m_exact@variance, m_capprox@variance, tolerance = 1e-5)
})

test_that("estimate() roxygen example runs under deriv_method = 'exact'", {
  # The shape of estimate()'s own runnable example, which is the ee_mean pattern
  # matrix(y - theta[1], nrow = 1) that previously errored under exact mode.
  psi <- function(theta) {
    y <- c(1, 2, 3, 4, 5)
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))

  m <- expect_no_error(estimate(m, deriv_method = "exact"))
  expect_equal(unname(m@theta), 3)
  # bread of the summed mean EE is 1 (d/dtheta of -sum(y - theta) is n; /n = 1)
  expect_equal(unname(m@bread), matrix(1))
})

test_that("estimate() fits a psi built with %*% under numerical methods", {
  # A user estimating function that forms the linear predictor with %*%
  # flows through the numerical bread path and recovers the OLS coefficients.
  set.seed(42)
  n <- 50
  X <- cbind(1, rnorm(n))
  y <- 1.5 + 2 * X[, 2] + rnorm(n)
  psi <- function(theta) {
    resid <- as.vector(y - X %*% theta)
    t(X * resid)
  }
  coefs <- unname(coef(lm(y ~ X[, 2])))

  se_capprox <- NULL
  for (method in c("capprox", "fapprox", "bapprox")) {
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0)),
      deriv_method = method
    )
    expect_equal(unname(m@theta), coefs, tolerance = 1e-6)
    expect_equal(dim(m@variance), c(2L, 2L))
    if (is.null(se_capprox)) {
      se_capprox <- sqrt(diag(m@variance))
    } else {
      expect_equal(sqrt(diag(m@variance)), se_capprox, tolerance = 1e-5)
    }
  }
})

# ---- the reduction a fit differentiates and solves ---------------------------
#
# A fit reduces the estimating functions to their row sums everywhere but the
# meat: once per solver evaluation, and once or twice per parameter while the
# bread is differentiated. Only the meat needs the per-observation contributions.
# An estimator whose `summed_equations` property holds those sums hands every
# reduction the closed form instead, and the meat takes the one full evaluation
# it always took.

reducer_fit_case <- function(n = 120) {
  set.seed(21)
  X <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  y <- as.vector(X %*% c(0.5, 1.5, -0.75)) + stats::rnorm(n)
  list(
    X = X,
    y = y,
    init = c(0, 0, 0),
    psi = function(theta) t(X * (y - as.vector(X %*% theta))),
    # The row sums of that psi as one matrix product. `t(X) %*% resid` carries
    # tangents through the registered methods, so the same reduction serves the
    # exact pass.
    summed = function(theta) {
      as.vector(t(X) %*% (y - as.vector(X %*% theta)))
    }
  )
}

counting_psi <- function(psi, counter) {
  function(theta) {
    counter$calls <- counter$calls + 1L
    psi(theta)
  }
}

test_that("a fit through a reduction evaluates psi only for the meat", {
  case <- reducer_fit_case()
  counter <- new.env(parent = emptyenv())

  counter$calls <- 0L
  MEstimator(
    stacked_equations = counting_psi(case$psi, counter),
    init = case$init,
    summed_equations = case$summed
  ) |>
    estimate()
  # The validation at the starting values, and the one the meat is built from.
  # The solver and the bread went through the reduction.
  expect_identical(counter$calls, 2L)

  counter$calls <- 0L
  MEstimator(
    stacked_equations = counting_psi(case$psi, counter),
    init = case$init
  ) |>
    estimate()
  # The same two, plus one per solver evaluation and two per parameter for the
  # central differences, which is what the property exists to avoid.
  expect_gt(counter$calls, 2L + 2L * length(case$init))
})

test_that("a fit through a reduction agrees with the fit through the matrix", {
  case <- reducer_fit_case()

  reduced <- MEstimator(
    stacked_equations = case$psi,
    init = case$init,
    summed_equations = case$summed
  ) |>
    estimate()
  matrix_route <- MEstimator(
    stacked_equations = case$psi,
    init = case$init
  ) |>
    estimate()

  expect_equal(coef(reduced), coef(matrix_route), tolerance = 1e-8)
  # The two reductions sum the same values in different orders, so their
  # difference quotients differ by rounding rather than exactly.
  expect_equal(vcov(reduced), vcov(matrix_route), tolerance = 1e-5)
})

test_that("a fit that supplies no reduction is the fit it always was", {
  case <- reducer_fit_case()

  expect_identical(
    MEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = NULL
    ) |>
      estimate() |>
      vcov(),
    MEstimator(stacked_equations = case$psi, init = case$init) |>
      estimate() |>
      vcov()
  )
})

test_that("a reduction of another system fails the fit", {
  case <- reducer_fit_case()
  # The second equation carries a term the estimating functions do not, so the
  # solver drives another system to zero and the bread would be the Jacobian of
  # that one while the meat is the cross-product of these.
  wrong <- function(theta) case$summed(theta) + c(0, 5, 0)

  err <- expect_error(
    MEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = wrong
    ) |>
      estimate(),
    class = "deli_summed_equations_disagree"
  )

  expect_match(conditionMessage(err), "summed_equations")
  expect_match(conditionMessage(err), "stacked_equations")
})

test_that("the fit's agreement check can be turned off", {
  case <- reducer_fit_case()
  wrong <- function(theta) case$summed(theta) + c(0, 5, 0)

  # What the check bought is then gone, and the fit reports the most it can
  # still see on its own: the estimating equations it was given are not solved
  # where the reduction sent it.
  expect_warning(
    fit <- MEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = wrong,
      check_summed_equations = FALSE
    ) |>
      estimate(),
    class = "deli_solver_not_converged"
  )

  expect_length(coef(fit), 3L)
})

test_that("a reduction serves a subset fit and an exact one", {
  case <- reducer_fit_case()

  subset_reduced <- MEstimator(
    stacked_equations = case$psi,
    init = case$init,
    subset = c(1L, 2L),
    summed_equations = case$summed
  ) |>
    estimate()
  subset_matrix <- MEstimator(
    stacked_equations = case$psi,
    init = case$init,
    subset = c(1L, 2L)
  ) |>
    estimate()
  expect_equal(coef(subset_reduced), coef(subset_matrix), tolerance = 1e-8)

  exact_reduced <- MEstimator(
    stacked_equations = case$psi,
    init = case$init,
    summed_equations = case$summed
  ) |>
    estimate(deriv_method = "exact")
  exact_matrix <- MEstimator(
    stacked_equations = case$psi,
    init = case$init
  ) |>
    estimate(deriv_method = "exact")
  expect_equal(vcov(exact_reduced), vcov(exact_matrix), tolerance = 1e-8)
})

test_that("a GMM fit takes the reduction for its objective and its bread", {
  case <- reducer_fit_case()
  counter <- new.env(parent = emptyenv())

  counter$calls <- 0L
  GMMEstimator(
    stacked_equations = counting_psi(case$psi, counter),
    init = case$init,
    summed_equations = case$summed
  ) |>
    estimate()
  # Just-identified, so no weight update runs: the validation at the starting
  # values and the one the meat is built from are the whole of it.
  expect_identical(counter$calls, 2L)

  counter$calls <- 0L
  GMMEstimator(
    stacked_equations = counting_psi(case$psi, counter),
    init = case$init
  ) |>
    estimate()
  expect_gt(counter$calls, 2L + 2L * length(case$init))
})

test_that("a GMM weight update keeps the per-observation contributions", {
  set.seed(31)
  counts <- stats::rpois(200, 3)
  psi <- function(theta) {
    rbind(counts - theta[1], (counts - theta[1])^2 - theta[1])
  }
  summed <- function(theta) {
    c(
      sum(counts) - 200 * theta[1],
      sum((counts - theta[1])^2) - 200 * theta[1]
    )
  }
  counter <- new.env(parent = emptyenv())
  counter$calls <- 0L

  # One weight-update pass, which a tolerance this wide settles after.
  GMMEstimator(
    stacked_equations = counting_psi(psi, counter),
    init = 1,
    summed_equations = summed,
    overid_maxiter = 1L,
    overid_tolerance = 1
  ) |>
    estimate()

  # The starting values, the one pass of the weight update, and the meat. The
  # weight matrix is the inverse of a cross-product, so that pass needs the
  # columns and the reduction cannot serve it.
  expect_identical(counter$calls, 3L)
})

test_that("a reduction that is not a function is refused at construction", {
  case <- reducer_fit_case()

  expect_error(
    MEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = c(1, 2, 3)
    ),
    class = "deli_summed_equations_error"
  )
  expect_error(
    GMMEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = c(1, 2, 3)
    ),
    class = "deli_summed_equations_error"
  )
})

test_that("the reduction properties default to no reduction and a check", {
  case <- reducer_fit_case()

  m <- MEstimator(stacked_equations = case$psi, init = case$init)
  expect_null(m@summed_equations)
  expect_true(m@check_summed_equations)

  g <- GMMEstimator(stacked_equations = case$psi, init = case$init)
  expect_null(g@summed_equations)
  expect_true(g@check_summed_equations)
})

test_that("a reduction of the wrong shape is refused before the solve", {
  case <- reducer_fit_case()
  # A reduction the bread cannot be differentiated from is refused wherever it
  # is supplied. Reaching the solver with one costs the diagnosis: the shape is
  # read by whatever base R arithmetic touches it first, several frames inside
  # a solver, against arguments the caller never named.
  bad <- list(
    short = function(theta) case$summed(theta)[1:2],
    text = function(theta) as.character(case$summed(theta)),
    list = function(theta) as.list(case$summed(theta))
  )

  for (label in names(bad)) {
    expect_error(
      MEstimator(
        stacked_equations = case$psi,
        init = case$init,
        summed_equations = bad[[label]]
      ) |>
        estimate(),
      class = "deli_summed_equations_error",
      info = label
    )
    expect_error(
      GMMEstimator(
        stacked_equations = case$psi,
        init = case$init,
        summed_equations = bad[[label]]
      ) |>
        estimate(),
      class = "deli_summed_equations_error",
      info = label
    )
  }
})

test_that("the refused shape names the starting values and both counts", {
  case <- reducer_fit_case()
  short <- function(theta) case$summed(theta)[1:2]

  err <- expect_error(
    MEstimator(
      stacked_equations = case$psi,
      init = case$init,
      summed_equations = short
    ) |>
      estimate(),
    class = "deli_summed_equations_error"
  )

  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_match(flat, "returned 2 values at the initial values", fixed = TRUE)
  expect_match(flat, "returns 3 estimating equations", fixed = TRUE)
  # A fit has no `theta` argument to point the caller at.
  expect_no_match(flat, "theta", fixed = TRUE)
})

test_that("the shape read at the starting values costs one reducer call", {
  case <- reducer_fit_case()
  counter <- new.env(parent = emptyenv())
  counter$psi <- 0L
  counter$summed <- 0L

  MEstimator(
    stacked_equations = function(theta) {
      counter$psi <- counter$psi + 1L
      case$psi(theta)
    },
    init = case$init,
    summed_equations = function(theta) {
      counter$summed <- counter$summed + 1L
      case$summed(theta)
    }
  ) |>
    estimate()

  # The shape read is a call to the reduction, not to the estimating functions,
  # so the count the validation and the meat account for is unchanged.
  expect_identical(counter$psi, 2L)
  expect_gt(counter$summed, 0L)
})

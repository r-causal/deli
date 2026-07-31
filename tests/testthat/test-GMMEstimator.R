# ---- GMMEstimator construction -----------------------------------------------

test_that("GMMEstimator constructs with valid inputs", {
  psi <- function(theta) {
    y <- c(1, 2, 3, 4, 5)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_s3_class(g, "deli::GMMEstimator")
  expect_equal(g@init, 0)
  expect_null(g@theta)
  expect_null(g@variance)
})

test_that("GMMEstimator stores init vector", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1))
  expect_equal(g@init, c(0, 1))
})

test_that("GMMEstimator accepts subset parameter", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = c(1L))
  expect_equal(g@subset, 1L)
})

test_that("GMMEstimator subset is sorted", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2], y * theta[3])
  }
  g <- GMMEstimator(
    stacked_equations = psi,
    init = c(0, 1, 1),
    subset = c(3L, 1L)
  )
  expect_equal(g@subset, c(1L, 3L))
})

test_that("GMMEstimator defaults are correct", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_null(g@subset)
  expect_null(g@finite_correction)
  expect_equal(g@overid_maxiter, 200L)
  expect_equal(g@overid_tolerance, 1e-9)
})

test_that("GMMEstimator accepts custom overid parameters", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(
    stacked_equations = psi,
    init = c(0),
    overid_maxiter = 20L,
    overid_tolerance = 1e-6
  )
  expect_equal(g@overid_maxiter, 20L)
  expect_equal(g@overid_tolerance, 1e-6)
})

test_that("GMMEstimator requires init to be numeric", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = "not numeric"),
    "numeric|double"
  )
})

test_that("GMMEstimator n_params is set from init length", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1))
  expect_equal(g@n_params, 2L)
})

# ---- GMMEstimator pre-estimation state ---------------------------------------

test_that("GMMEstimator theta is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_null(g@theta)
})

test_that("GMMEstimator variance is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_null(g@variance)
})

test_that("GMMEstimator weight_matrix is NULL before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_null(g@weight_matrix)
})

# ---- GMMEstimator constructor validation -------------------------------------
#
# These tests pin construction-time validation of finite_correction, init,
# subset, and the GMM-specific overid_maxiter and overid_tolerance arguments.
# finite_correction, init, and subset follow the same contract as MEstimator:
# finite_correction is NULL or the HC1 string, init must contain at least one
# value, and subset must be whole-number parameter indices between 1 and the
# number of parameters. overid_maxiter and overid_tolerance must each be a
# single number. Note on parity: the estimator runs zero over-identification
# iterations when overid_maxiter is below 1 and never converges early when
# overid_tolerance is non-positive, matching Python delicatessen, so the
# validators reject non-numeric, missing, and non-scalar values rather than
# non-positive ones.

test_that("GMMEstimator rejects an unsupported finite_correction string", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = "nonsense"
    ),
    "HC1"
  )
})

test_that("GMMEstimator rejects a non-character finite_correction", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = 5
    ),
    "HC1"
  )
})

test_that("GMMEstimator accepts the supported finite_correction values", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_no_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = "HC1"
    )
  )
  expect_no_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      finite_correction = NULL
    )
  )
})

test_that("GMMEstimator rejects zero-length init", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = numeric(0)),
    "at least one"
  )
})

test_that("GMMEstimator rejects a subset index outside the parameter range", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = 5L),
    "subset"
  )
})

test_that("GMMEstimator rejects a zero subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = 0L),
    "subset"
  )
})

test_that("GMMEstimator rejects a negative subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = -1L),
    "subset"
  )
})

test_that("GMMEstimator rejects a non-numeric subset", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = "a"),
    "subset"
  )
})

test_that("GMMEstimator rejects a fractional subset index", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = 1.5),
    "subset"
  )
})

test_that("GMMEstimator accepts a valid in-range subset", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_no_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = c(1L, 2L))
  )
})

test_that("GMMEstimator rejects a subset with duplicated indices", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0, 1), subset = c(2L, 2L)),
    "duplicated"
  )
})

test_that("GMMEstimator rejects a non-numeric overid_maxiter", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  # Target message: overid_maxiter must be a single number.
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0), overid_maxiter = "x"),
    "single"
  )
})

test_that("GMMEstimator rejects a non-scalar overid_maxiter", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      overid_maxiter = c(1L, 2L)
    ),
    "single"
  )
})

test_that("GMMEstimator rejects a non-numeric overid_tolerance", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  # Target message: overid_tolerance must be a single number.
  expect_error(
    GMMEstimator(stacked_equations = psi, init = c(0), overid_tolerance = "x"),
    "single"
  )
})

test_that("GMMEstimator rejects a non-scalar overid_tolerance", {
  psi <- function(theta) {
    matrix(c(1, 2, 3) - theta[1], nrow = 1)
  }
  expect_error(
    GMMEstimator(
      stacked_equations = psi,
      init = c(0),
      overid_tolerance = c(1e-9, 1e-8)
    ),
    "single"
  )
})

# ---- estimate() just-identified (same as MEstimator) -------------------------

test_that("GMMEstimator finds correct theta for mean EE", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  # GMM should give same results as M-estimator for just-identified
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(g@theta, m@theta, tolerance = 1e-4)
})

test_that("GMMEstimator finds correct variance for mean EE", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(g@variance[1, 1], m@variance[1, 1], tolerance = 1e-4)
})

test_that("GMMEstimator finds correct theta for mean+variance EE", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(g@theta, m@theta, tolerance = 1e-4)
})

test_that("GMMEstimator computes correct variance for mean+variance EE", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m <- estimate(m)

  expect_equal(diag(g@variance), diag(m@variance), tolerance = 1e-4)
})

# ---- estimate() sets properties correctly ------------------------------------

test_that("estimate() sets n_obs for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  expect_equal(g@n_obs, 5L)
})

test_that("estimate() populates bread and meat for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  expect_true(!is.null(g@bread))
  expect_true(!is.null(g@meat))
  expect_true(is.matrix(g@bread))
  expect_true(is.matrix(g@meat))
})

test_that("estimate() sets weight_matrix for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  expect_true(!is.null(g@weight_matrix))
  # For just-identified, weight_matrix should be identity
  expect_equal(g@weight_matrix, diag(1))
})

test_that("estimate() returns a GMMEstimator object", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  result <- estimate(g)

  expect_s3_class(result, "deli::GMMEstimator")
})

# ---- Vector-return psi (single parameter) ------------------------------------
#
# A psi returning a plain vector (no dim attribute) is the valid single-parameter
# form and already works in MEstimator, which reshapes it to a 1-by-n matrix
# before forming the meat. The GMM path must apply the same reshape so the meat
# cross-product is 1-by-1 rather than n-by-n.

test_that("GMMEstimator handles a vector-return psi like MEstimator", {
  y <- c(1, 2, 3, 4, 5, 6)
  psi <- function(theta) y - theta[1]

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m)

  expect_equal(unname(g@theta), mean(y), tolerance = 1e-5)
  expect_equal(unname(g@theta), unname(m@theta), tolerance = 1e-5)
  expect_equal(
    unname(g@variance[1, 1]),
    unname(m@variance[1, 1]),
    tolerance = 1e-5
  )
})

# ---- Over-identified estimation ----------------------------------------------

test_that("GMMEstimator handles over-identified equations", {
  # 1 parameter, 2 equations (over-identified)
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    rbind(
      y - theta[1], # mean equation
      (y - theta[1])^2 - 2 # variance equation (moment condition)
    )
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(2))
  g <- estimate(g)

  # Should find a reasonable estimate near the mean
  expect_equal(unname(g@theta[1]), mean(y), tolerance = 0.5)
  expect_true(!is.null(g@variance))

  # Weight matrix should NOT be identity for over-identified
  expect_false(isTRUE(all.equal(g@weight_matrix, diag(2))))
})

test_that("GMMEstimator errors on under-identified equations", {
  # 2 parameters, 1 equation (under-identified)
  y <- c(1, 2, 3, 4, 5)

  psi <- function(theta) {
    matrix(y - theta[1] - theta[2], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 0))
  expect_error(estimate(g), "less than or equal")
})

test_that("GMMEstimator reports under-identification ahead of non-finite moments", {
  # 3 parameters against 2 equations, and non-finite at the starting values
  # because log(theta[2]) is -Inf there. No starting values make an
  # under-identified system solvable, so the shortfall is what gets reported.
  y <- c(1, 2, 3, 4, 5)

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - log(theta[2])
    )
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 0, 0))
  err <- expect_error(estimate(g), class = "deli_psi_shape_error")
  expect_false(grepl("non-finite", conditionMessage(err), fixed = TRUE))
})

# The two solver diagnostics below are cli warnings (rlang_warning), matching the
# other convergence conditions in estimate(). The message substrings are the
# stable part callers and tests key on.

test_that("GMMEstimator signals over-identification non-convergence as a cli warning", {
  # Over-identified system capped at a single updating step with a tolerance it
  # cannot meet, so the iterative GMM loop terminates without converging.
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - 2
    )
  }

  g <- GMMEstimator(
    stacked_equations = psi,
    init = c(2),
    overid_maxiter = 1L,
    overid_tolerance = 1e-12
  )

  cnd <- rlang::catch_cnd(estimate(g), classes = "warning")
  expect_s3_class(cnd, "rlang_warning")
  expect_match(conditionMessage(cnd), "iterative GMM updating")
})

# ---- the default over-identification budget ----------------------------------
#
# The two-step weight update converges linearly, not quadratically. Measured on
# well-specified linear IV systems the per-iteration error ratio runs between
# 0.42 and 0.88, so reaching the default `overid_tolerance` of 1e-9 takes tens
# of updates rather than the handful a budget of ten allows. Bisecting the
# budget over a hundred seeds of the construction below puts the median at 17
# updates and the worst case at 42, and four fifths of those seeds need more
# than ten. A well-specified, well-identified fit that warns is reporting the
# budget rather than anything about the data, so the default has to clear the
# range the data actually needs.
#
# Three instruments for an intercept and a treatment effect: three moment
# conditions for two parameters, over-identified by one. The instruments are
# strong and the outcome model is correct, so nothing here is a specification
# failure and the budget is the only thing that can stop the fit short.

make_overid_iv_psi <- function(seed, n = 500) {
  set.seed(seed)
  instruments <- cbind(
    stats::rbinom(n, 1, 0.5),
    stats::rnorm(n),
    stats::rnorm(n)
  )
  confounder <- stats::rnorm(n)
  treatment <- as.numeric(instruments %*% rep(0.4, 3)) +
    confounder +
    stats::rnorm(n)
  outcome <- 1 + 2 * treatment - confounder + stats::rnorm(n)
  design <- cbind(1, treatment)
  function(theta) {
    t(instruments * as.numeric(outcome - design %*% theta))
  }
}

overid_iv_init <- c(intercept = 0, effect = 0)

# Seeds of that construction whose updating loop needs 32, 26 and 24 updates to
# settle, measured by bisecting `overid_maxiter` against the real fit. All three
# are inside any budget the measurements support and outside a budget of ten.
slow_overid_seeds <- c(3L, 11L, 12L)

test_that("the default overid_maxiter clears what a well-specified fit needs", {
  # The measured worst case over a hundred seeds is 42 updates, and a weakly
  # identified system can want several hundred, so the default is set well above
  # the range a well-identified fit occupies rather than at the edge of it.
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  expect_identical(
    GMMEstimator(stacked_equations = psi, init = c(0))@overid_maxiter,
    200L
  )
  expect_identical(formals(gmm_estimate.default)$overid_maxiter, quote(200L))
  expect_identical(formals(gmm_estimate.formula)$overid_maxiter, quote(200L))
})

test_that("a well-specified over-identified fit is silent at the default budget", {
  for (seed in slow_overid_seeds) {
    psi <- make_overid_iv_psi(seed)
    expect_no_warning(
      gmm_estimate(stacked_equations = psi, init = overid_iv_init)
    )
  }
})

test_that("the default budget reaches the same fit as a generous one", {
  # Silence on its own would be satisfied by loosening the tolerance, so the
  # point the default lands on is compared with the point a budget nobody
  # doubts lands on, to the last bit rather than to a tolerance.
  for (seed in slow_overid_seeds) {
    psi <- make_overid_iv_psi(seed)
    at_default <- expect_no_warning(
      gmm_estimate(stacked_equations = psi, init = overid_iv_init)
    )
    generous <- gmm_estimate(
      stacked_equations = psi,
      init = overid_iv_init,
      overid_maxiter = 400L
    )
    expect_identical(at_default@theta, generous@theta)
    expect_identical(at_default@variance, generous@variance)
    expect_identical(at_default@weight_matrix, generous@weight_matrix)
  }
})

test_that("the over-identified IV system of the examples is silent at the default", {
  # The system the `gmm_estimate()` examples fit, which needs 26 updates and so
  # passed an explicit `overid_maxiter` to stay quiet.
  set.seed(42)
  n <- 200
  z1 <- stats::rbinom(n, 1, 0.5)
  z2 <- stats::rnorm(n)
  confounder <- stats::rnorm(n)
  treatment <- 0.5 * z1 + 0.3 * z2 + confounder + stats::rnorm(n)
  outcome <- 2 * treatment - confounder + stats::rnorm(n)
  psi <- function(theta) {
    residual <- outcome - theta[1] * treatment
    rbind(z1 * residual, z2 * residual)
  }
  expect_no_warning(gmm_estimate(
    stacked_equations = psi,
    init = c(effect = 0)
  ))
})

test_that("an updating loop that does not settle still warns at the default", {
  # Raising the budget cannot cure weak identification: the same construction
  # with instruments a tenth as strong takes several hundred updates to settle,
  # or never settles, and saying so is the warning's job.
  set.seed(1)
  n <- 300
  instruments <- cbind(
    stats::rbinom(n, 1, 0.5),
    stats::rnorm(n),
    stats::rnorm(n)
  )
  confounder <- stats::rnorm(n)
  treatment <- as.numeric(instruments %*% rep(0.05, 3)) +
    confounder +
    stats::rnorm(n)
  outcome <- 1 + 2 * treatment - confounder + stats::rnorm(n)
  design <- cbind(1, treatment)
  psi <- function(theta) {
    t(instruments * as.numeric(outcome - design %*% theta))
  }
  cnd <- rlang::catch_cnd(
    gmm_estimate(stacked_equations = psi, init = overid_iv_init),
    classes = "warning"
  )
  expect_s3_class(cnd, "rlang_warning")
  expect_match(conditionMessage(cnd), "iterative GMM updating")
})

test_that("GMMEstimator signals solver non-convergence as a cli warning", {
  # A single optimizer iteration cannot reach the mean, so optim reports a
  # non-zero convergence code.
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))

  cnd <- rlang::catch_cnd(estimate(g, maxiter = 1), classes = "warning")
  expect_s3_class(cnd, "rlang_warning")
  expect_match(conditionMessage(cnd), "did not converge")
})

# ---- Inference methods work for GMMEstimator ---------------------------------

test_that("confidence_intervals works for GMMEstimator", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  ci <- confidence_intervals(g)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  expect_equal(colnames(ci), c("lower", "upper"))
  # Mean should be inside interval
  expect_true(ci[1, "lower"] < mean(y))
  expect_true(ci[1, "upper"] > mean(y))
})

test_that("z_scores works for GMMEstimator", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  z <- z_scores(g)
  expect_true(is.numeric(z))
  expect_equal(length(z), 1)
})

test_that("p_values works for GMMEstimator", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  p <- p_values(g)
  expect_true(is.numeric(p))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("s_values works for GMMEstimator", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  s <- s_values(g)
  expect_true(is.numeric(s))
  expect_true(all(s >= 0))
})

# ---- Print and summary work for GMMEstimator ---------------------------------

test_that("print works for GMMEstimator before estimation", {
  psi <- function(theta) {
    y <- c(1, 2, 3)
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  expect_snapshot(print(g))
})

test_that("print works for GMMEstimator after estimation", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)
  expect_snapshot(print(g))
})

test_that("summary works for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  g <- estimate(g)

  s <- summary(g)
  expect_s3_class(s, "deli::EstimatorSummary")
  expect_equal(s@theta, g@theta)
})

# ---- Custom solver for GMMEstimator ------------------------------------------

test_that("GMMEstimator accepts custom solver function", {
  y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  custom_solver <- function(stacked_equations, init) {
    result <- stats::optim(
      par = init,
      fn = stacked_equations,
      method = "Nelder-Mead"
    )
    result$par
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 1))
  g <- estimate(g, solver = custom_solver)

  expect_equal(unname(g@theta[1]), mean(y), tolerance = 1e-3)
})

# ---- Regression with GMMEstimator -------------------------------------------

test_that("GMMEstimator works for linear regression EE", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 3 * x + rnorm(n)
  X <- cbind(1, x)

  psi <- function(theta) {
    resid <- y - X %*% theta
    # Each row is an estimating equation, each column is an observation
    t(X) * matrix(resid, nrow = ncol(X), ncol = n, byrow = TRUE)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 0))
  g <- estimate(g)

  m <- MEstimator(stacked_equations = psi, init = c(0, 0))
  m <- estimate(m)

  # GMM and M-estimator should give similar results for just-identified
  expect_equal(g@theta, m@theta, tolerance = 1e-3)
  expect_equal(diag(g@variance), diag(m@variance), tolerance = 1e-3)
})

# ---- what the bread of a GMM fit says about identification -------------------
#
# A GMM fit judges its moments, through the J-statistic where it is
# over-identified and through the readings of the returned point where it is
# not, but nothing asked the bread whether the parameters can be told apart. A
# rank-deficient just-identified system therefore came back silent, where the
# same system fitted as an M-estimator reports that its parameters are not
# identified. The reading is the M path's, reused.

rank_deficient_gmm_data <- function() {
  set.seed(11)
  n <- 60
  x <- stats::rnorm(n)
  # The third column repeats the second, so no value of the parameters is
  # distinguishable from another that trades between the two.
  list(X = cbind(1, x, x), y = 1 + 2 * x + stats::rnorm(n))
}

test_that("a rank-deficient just-identified GMM fit reports its parameters", {
  d <- rank_deficient_gmm_data()
  psi <- function(theta) t(d$X * as.vector(d$y - d$X %*% theta))

  g <- GMMEstimator(stacked_equations = psi, init = c(0, 0, 0))
  cnd <- rlang::catch_cnd(estimate(g), classes = "warning")

  expect_s3_class(cnd, "deli_solver_not_converged")
  expect_match(conditionMessage(cnd), "not identified")
})

test_that("a healthy over-identified GMM fit stays silent about identification", {
  # The bread of an over-identified fit is rectangular by design, so the
  # question its rank answers is whether the columns are independent rather
  # than whether the matrix is square. Reading squareness would report every
  # over-identified fit.
  set.seed(12)
  n <- 200
  z1 <- stats::rbinom(n, 1, 0.5)
  z2 <- stats::rnorm(n)
  u <- stats::rnorm(n)
  a <- 0.5 * z1 + 0.3 * z2 + u + stats::rnorm(n)
  y <- 2 * a - u + stats::rnorm(n)

  psi <- function(theta) {
    resid <- y - theta[1] * a
    rbind(z1 * resid, z2 * resid)
  }

  g <- GMMEstimator(stacked_equations = psi, init = 0)
  expect_no_warning(estimate(g))
  expect_equal(unname(coef(estimate(g))), 2, tolerance = 0.3)
})

# ---- the moment covariance the weight update cannot invert -------------------

test_that("a singular moment covariance is classed under allow_pinv = FALSE", {
  # The update takes the inverse of the moment covariance, and under
  # `allow_pinv = FALSE` it took it with a bare `solve()`, so a covariance with
  # no inverse surfaced as LAPACK's "system is exactly singular" against an
  # argument no caller wrote. The bread has carried a class for this since
  # `check_bread_invertible()`; the meat had none.
  set.seed(13)
  n <- 120
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  a <- z1 + 0.5 * z2 + stats::rnorm(n)
  y <- 2 * a + stats::rnorm(n)

  # The second moment condition repeats the first exactly, so the covariance of
  # the three is singular.
  psi <- function(theta) {
    resid <- y - theta[1] * a
    rbind(z1 * resid, z1 * resid, z2 * resid)
  }

  g <- GMMEstimator(stacked_equations = psi, init = 0)
  err <- expect_error(
    estimate(g, allow_pinv = FALSE),
    class = "deli_meat_not_invertible"
  )
  flat <- gsub("\\s+", " ", conditionMessage(err))
  expect_match(flat, "allow_pinv", fixed = TRUE)
  expect_match(flat, "rank is 2 of 3", fixed = TRUE)
  # LAPACK's own account of the failure named an argument no caller wrote.
  expect_false(grepl("reciprocal condition number", flat, fixed = TRUE))
})

test_that("the dependent-moment report degrades where qr sees full rank", {
  # `warn_dependent_moments()` is reached whenever `solve()` fails, and it
  # assumed that a failed solve means a rank-deficient factorization. A
  # covariance whose columns are orthogonal but whose scales differ by more than
  # the reciprocal condition tolerance defeats `solve()` while `qr()` still
  # reports full rank, which left the report naming a rank of k out of k and
  # listing no condition at all.
  meat <- diag(c(1, 1e-20))
  expect_identical(qr(meat)$rank, 2L)

  cnd <- rlang::catch_cnd(warn_dependent_moments(meat), classes = "warning")
  expect_s3_class(cnd, "deli_gmm_moments_dependent")
  flat <- gsub("\\s+", " ", conditionMessage(cnd))
  expect_match(flat, "poorly conditioned", fixed = TRUE)
  expect_match(flat, "full rank", fixed = TRUE)
  # The sentence that names the dependent conditions has none to name here.
  expect_false(grepl("accounted for by the others, so", flat, fixed = TRUE))
})

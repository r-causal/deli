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
  expect_equal(g@overid_maxiter, 10L)
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

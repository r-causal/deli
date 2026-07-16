# ---- Helper ------------------------------------------------------------------

make_fitted_mean <- function(solver = "rootSolve") {
  ref <- load_fixture("ee_mean")
  y <- ref$y
  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }
  m <- MEstimator(stacked_equations = psi, init = c(0))
  estimate(m, solver = solver)
}

# ---- Default rootSolve solver ------------------------------------------------

test_that("rootSolve solver returns correct theta and variance", {
  ref <- load_fixture("ee_mean")
  m <- make_fitted_mean(solver = "rootSolve")

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
  expect_equal(unname(m@variance[1, 1]), ref$variance[[1]], tolerance = 1e-4)
})

# ---- nleqslv solver ----------------------------------------------------------

test_that("nleqslv solver gives same results as rootSolve for ee_mean", {
  ref <- load_fixture("ee_mean")
  m_nleqslv <- make_fitted_mean(solver = "nleqslv")

  expect_equal(unname(m_nleqslv@theta), ref$theta, tolerance = 1e-5)
  expect_equal(unname(m_nleqslv@variance[1, 1]), ref$variance[[1]], tolerance = 1e-4)
})

# ---- Custom solver function --------------------------------------------------

test_that("custom solver function works", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  # Simple custom solver that uses stats::optim to minimize sum of squares

  my_solver <- function(stacked_equations, init) {
    obj <- function(theta) {
      sum(stacked_equations(theta)^2)
    }
    result <- optim(init, obj, method = "BFGS")
    result$par
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  m <- estimate(m, solver = my_solver)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
})

# ---- Invalid solver inputs ---------------------------------------------------

test_that("invalid solver string errors", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m, solver = "nonexistent"))
})

test_that("invalid solver type (non-string, non-function) errors", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(estimate(m, solver = 42))
  expect_error(estimate(m, solver = TRUE))
  expect_error(estimate(m, solver = list()))
})

test_that("unknown solver string reports it is not supported", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(
    estimate(m, solver = "definitely_not_a_solver"),
    regexp = "not supported"
  )
})

# ---- Levenberg-Marquardt solver (solver = "lm") ------------------------------
#
# These tests pin the API surface and numerical behaviour of the
# Levenberg-Marquardt solver. In Python delicatessen it is the default,
# reached through scipy.optimize.root(method = "lm"); the deli surface
# exposes it as solver = "lm". The tracked implementation is bd-rdyz.
#
# Fixtures ee_solver_lm_linear and ee_solver_lm_logistic come from Python
# solver = "lm" (fixtures/generate_solver_lm_fixtures.py). The linear
# problem is well-conditioned so every root-finder agrees to tight
# tolerance; the logistic problem is a nonlinear system that exercises the
# least-squares behaviour of Levenberg-Marquardt.

# Build the linear regression estimating equation from a fixture.
fit_linear <- function(fixture, solver) {
  ref <- load_fixture(fixture)
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    residuals <- y - X %*% theta
    t(X * as.numeric(residuals))
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  estimate(m, solver = solver)
}

# Build the logistic regression estimating equation from a fixture.
fit_logistic <- function(fixture, solver) {
  ref <- load_fixture(fixture)
  X <- ref$X
  y <- ref$y

  inverse_logit <- function(x) 1 / (1 + exp(-x))

  psi <- function(theta) {
    residuals <- y - inverse_logit(X %*% theta)
    t(X * as.numeric(residuals))
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  estimate(m, solver = solver)
}

# Standard errors from a fitted estimator, for comparisons on the SE scale.
param_se <- function(m) sqrt(diag(m@variance))

test_that("lm solver matches Python solver='lm' for linear regression", {
  m <- fit_linear("ee_solver_lm_linear", solver = "lm")
  expect_python_match(m, "ee_solver_lm_linear")
})

test_that("lm solver matches Python solver='lm' for logistic regression", {
  m <- fit_logistic("ee_solver_lm_logistic", solver = "lm")
  expect_python_match(m, "ee_solver_lm_logistic")
})

test_that("lm solver agrees with rootSolve on the linear problem", {
  m_lm <- fit_linear("ee_solver_lm_linear", solver = "lm")
  m_rs <- fit_linear("ee_solver_lm_linear", solver = "rootSolve")

  expect_equal(unname(m_lm@theta), unname(m_rs@theta), tolerance = 1e-6)
  expect_equal(unname(param_se(m_lm)), unname(param_se(m_rs)), tolerance = 1e-6)
})

test_that("lm solver agrees with rootSolve on the logistic problem", {
  m_lm <- fit_logistic("ee_solver_lm_logistic", solver = "lm")
  m_rs <- fit_logistic("ee_solver_lm_logistic", solver = "rootSolve")

  expect_equal(unname(m_lm@theta), unname(m_rs@theta), tolerance = 1e-6)
  expect_equal(unname(param_se(m_lm)), unname(param_se(m_rs)), tolerance = 1e-6)
})

# The existing rootSolve solver already reproduces the Python solver='lm'
# roots on these problems, which confirms the fixtures independently of the
# lm backend. Comparisons on theta use a relative tolerance, and variance is
# checked on the SE scale where all magnitudes are away from zero.
test_that("rootSolve reproduces the Python solver='lm' linear fixture", {
  ref <- load_fixture("ee_solver_lm_linear")
  m <- fit_linear("ee_solver_lm_linear", solver = "rootSolve")

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
  expect_equal(
    unname(param_se(m)),
    sqrt(diag(as.matrix(ref$variance))),
    tolerance = 1e-5
  )
})

test_that("rootSolve reproduces the Python solver='lm' logistic fixture", {
  ref <- load_fixture("ee_solver_lm_logistic")
  m <- fit_logistic("ee_solver_lm_logistic", solver = "rootSolve")

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-5)
  expect_equal(
    unname(param_se(m)),
    sqrt(diag(as.matrix(ref$variance))),
    tolerance = 1e-5
  )
})

# ---- rootSolve non-convergence detection -------------------------------------
#
# rootSolve::multiroot reports success even when the estimating-function
# Jacobian is singular at the starting values: it returns those values
# unchanged with a non-zero score and no warning. The lm and nleqslv branches
# already warn when they fail to converge; the default rootSolve branch now
# surfaces the same signal so the silent failure is visible. The mroz 2SLS
# problem started from zero is the reference reproducer: Python delicatessen
# defaults to Levenberg-Marquardt, which escapes the singular start, so a user
# who matches the Python init otherwise receives a silently wrong answer.

test_that("rootSolve warns when the Jacobian is singular at the start", {
  d <- mroz[stats::complete.cases(mroz), ]
  d$intercept <- 1
  psi <- function(theta) {
    ee_2sls(
      theta,
      y = d$lwage,
      A = d$educ,
      Z = as.matrix(d[, c("motheduc", "fatheduc", "huseduc")]),
      W = as.matrix(d[, c("intercept", "exper", "expersq")])
    )
  }

  expect_warning(
    m <- m_estimate(psi, init = rep(0, 10)),
    "rootSolve did not converge"
  )
  # The returned point is still the spurious root rootSolve reported; the fix
  # only surfaces the warning, it does not change the returned values.
  expect_equal(unname(m@theta), rep(0, 10))
})

test_that("a well-conditioned rootSolve fit does not warn", {
  m <- fit_linear("ee_solver_lm_linear", solver = "rootSolve")
  expect_no_warning(fit_linear("ee_solver_lm_linear", solver = "rootSolve"))
  expect_equal(
    unname(m@theta),
    load_fixture("ee_solver_lm_linear")$theta,
    tolerance = 1e-5
  )
})

# The non-smooth median equation in ee_positive_mean_deviation legitimately
# settles at a point where the summed score is not zero, yet it is not a
# spurious singular-Jacobian return. Initializing at the Python solution, the
# quality check must stay quiet there so the convergence warning does not fire
# on a genuine (if non-differentiable) fit.
test_that("rootSolve stays quiet on the non-smooth positive mean deviation", {
  ref <- load_fixture("ee_positive_mean_deviation")
  y <- ref$y
  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }
  expect_no_warning({
    m <- MEstimator(stacked_equations = psi, init = ref$theta)
    m <- estimate(m)
  })
})

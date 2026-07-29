# ---- Helper ------------------------------------------------------------------

# The reading of a returned point that uses nothing but the contributions of the
# equations there, which is one of the two readings unsolved_point() combines.
# Naming it separately is what lets the tests below pin each reading on its own.
is_root <- function(ef, theta) {
  is.null(unsolved_equation(ef, theta))
}

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
  expect_equal(
    unname(m_nleqslv@variance[1, 1]),
    ref$variance[[1]],
    tolerance = 1e-4
  )
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

# The body of the MEstimator method sits in an internal worker, which is a frame
# no user can name and which appears in no man page. The error therefore has to
# report the method's own frame instead of the worker's.
test_that("invalid solver type reports a frame the user can name", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  m <- MEstimator(stacked_equations = psi, init = c(0))
  err <- tryCatch(estimate(m, solver = 42), error = function(e) e)
  reported <- paste(deparse(conditionCall(err)), collapse = " ")

  expect_match(reported, "method(estimate, deli::MEstimator)", fixed = TRUE)
  expect_false(grepl("estimate_m_estimator", reported, fixed = TRUE))
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
# exposes it as solver = "lm".
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
  skip_if_not_installed("minpack.lm")
  m <- fit_linear("ee_solver_lm_linear", solver = "lm")
  expect_python_match(m, "ee_solver_lm_linear")
})

test_that("lm solver matches Python solver='lm' for logistic regression", {
  skip_if_not_installed("minpack.lm")
  m <- fit_logistic("ee_solver_lm_logistic", solver = "lm")
  expect_python_match(m, "ee_solver_lm_logistic")
})

test_that("lm solver agrees with rootSolve on the linear problem", {
  skip_if_not_installed("minpack.lm")
  m_lm <- fit_linear("ee_solver_lm_linear", solver = "lm")
  m_rs <- fit_linear("ee_solver_lm_linear", solver = "rootSolve")

  expect_equal(unname(m_lm@theta), unname(m_rs@theta), tolerance = 1e-6)
  expect_equal(unname(param_se(m_lm)), unname(param_se(m_rs)), tolerance = 1e-6)
})

test_that("lm solver agrees with rootSolve on the logistic problem", {
  skip_if_not_installed("minpack.lm")
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

# A capped iteration budget is a distinct non-convergence mode: rootSolve stops
# at a partially-solved point where the per-equation contributions still cancel
# well, so the singular-Jacobian cancellation heuristic stays quiet. rootSolve
# itself emits a "steady-state not reached" warning that the suppression around
# the call discards, so the maxiter exhaustion must be surfaced separately. The
# lm and nleqslv branches already warn on their own non-convergence codes.
test_that("rootSolve warns when the iteration budget is exhausted", {
  ref <- load_fixture("ee_solver_lm_logistic")
  X <- ref$X
  y <- ref$y
  inverse_logit <- function(x) 1 / (1 + exp(-x))
  psi <- function(theta) {
    residuals <- y - inverse_logit(X %*% theta)
    t(X * as.numeric(residuals))
  }
  fit <- MEstimator(stacked_equations = psi, init = ref$init)

  expect_warning(
    estimate(fit, maxiter = 2),
    "did not converge"
  )
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

# ---- rootSolve divergence detection ------------------------------------------
#
# multiroot's own convergence test can fail while the returned list still looks
# like a success: it carries no status flag, and the "steady-state not reached"
# warning rootSolve raises alongside it was discarded by the suppression around
# the call. A beta regression started from a scale parameter that is far too
# small is the reference reproducer. multiroot walks the location parameters out
# to 1e7 and stops with a largest score of 36 against a tolerance of 1e-9.
# Neither existing check can see it: the iteration budget is nowhere near
# exhausted (six iterations of five thousand), and 17% of the per-observation
# contributions still cancel, which keeps the cancellation heuristic quiet.
# Nor does comparing scores help, because the mean absolute score falls from
# 36.9 at the starting values to 24.3 at the returned point, so by that measure
# the solver made progress while sending the estimates to nonsense.

beta_runaway_psi <- function() {
  set.seed(42)
  n <- 50
  w <- stats::rnorm(n)
  mu <- stats::plogis(0.5 + 0.3 * w)
  y <- stats::rbeta(n, mu * 10, (1 - mu) * 10)
  x <- cbind(1, w)
  function(theta) ee_beta_regression(theta, X = x, y = y)
}

# Every warning raised while evaluating `expr`, with each one muffled so it
# cannot escape into the test summary. Used where a single call raises more than
# one warning, or where the count itself is the assertion.
collect_warnings <- function(expr) {
  seen <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  seen
}

test_that("rootSolve warns when it walks away from the starting values without converging", {
  psi <- beta_runaway_psi()
  expect_warning(
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0, log(10)))),
    class = "deli_solver_not_converged"
  )
  # Surfacing the failure does not change what the solver returned.
  expect_gt(max(abs(m@theta)), 1e6)
})

test_that("a diverging rootSolve solve raises exactly one warning", {
  psi <- beta_runaway_psi()
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = c(0, 0, log(10))))
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
})

test_that("a diverging rootSolve solve whose point is also unsolved warns once", {
  set.seed(42)
  n <- 50
  w <- stats::rnorm(n)
  mu <- stats::plogis(0.5 + 0.3 * w)
  y <- stats::rbeta(n, mu * 10, (1 - mu) * 10)
  x <- cbind(1, w)
  # The runaway beta stack with a fourth equation no parameter value solves, so
  # the returned point fails the readings as well as leaving multiroot
  # reporting its convergence test failed. Both paths would report this solve
  # and only one of them may.
  psi <- function(theta) {
    rbind(
      ee_beta_regression(theta[1:3], X = x, y = y),
      rep(theta[4]^2 + 1, n)
    )
  }
  seen <- collect_warnings({
    m <- estimate(MEstimator(
      stacked_equations = psi,
      init = c(0, 0, log(10), 1)
    ))
  })
  expect_length(seen, 1L)
  expect_match(
    conditionMessage(seen[[1]]),
    "did not converge to a root",
    fixed = TRUE
  )
  theta <- unname(m@theta)
  bread <- compute_bread(psi, theta, "capprox", 1e-9) / n
  expect_false(is.null(unsolved_point(psi(theta), theta, bread)))
})

test_that("an exhausted rootSolve budget warns once and is not judged again", {
  set.seed(3)
  y1 <- stats::rnorm(200, mean = 4)
  y2 <- stats::rnorm(200, mean = 2)
  psi <- function(theta) {
    rbind(
      y1 - theta[1],
      y2 - theta[2],
      rep(theta[1] / theta[2] - theta[3], 200)
    )
  }
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(1, 1, 1)),
      maxiter = 1
    )
  })
  expect_length(seen, 1L)
  expect_match(conditionMessage(seen[[1]]), "within 1 iteration", fixed = TRUE)
  # One iteration leaves the ratio row far from solved, so the returned point
  # would be reported a second time if the exhausted budget did not already
  # count as a warning about this solve.
  theta <- unname(m@theta)
  bread <- compute_bread(psi, theta, "capprox", 1e-9) / 200
  expect_false(is.null(unsolved_point(psi(theta), theta, bread)))
})

test_that("the same beta regression stays quiet from starting values that converge", {
  psi <- beta_runaway_psi()
  expect_no_warning({
    m <- estimate(MEstimator(
      stacked_equations = psi,
      init = c(0.5, 0.3, log(10))
    ))
  })
  expect_equal(unname(m@theta), c(0.46209, 0.278302, 2.55663), tolerance = 1e-4)
})

# A score of NaN reaches the same check. Comparing it against the score floor
# yields NA, which used to abort the whole fit with "missing value where
# TRUE/FALSE needed" rather than reporting a failed solve.
test_that("rootSolve reports a non-finite score instead of aborting", {
  psi <- beta_runaway_psi()
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = c(1, 1, 2)))
  )
  expect_true(
    any(vapply(seen, inherits, logical(1), "deli_solver_not_converged"))
  )
})

# ---- fits that must stay quiet -----------------------------------------------
#
# Judging a returned root by the size of its residual cannot separate a failed
# solve from a correct one: these five fits are all correct, all currently
# quiet, and all sit above any residual threshold that catches the divergence
# above. They are locked here so that no later convergence check reintroduces a
# false alarm on them. The positive-mean-deviation fixture is locked by
# "rootSolve stays quiet on the non-smooth positive mean deviation" above.

test_that("ee_positive_mean_deviation at n = 200 does not warn", {
  set.seed(2024)
  y <- 0.6 + stats::rnorm(200)
  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }
  expect_no_warning(
    estimate(MEstimator(stacked_equations = psi, init = c(1, stats::median(y))))
  )
})

test_that("ee_percentile at the median of a tied sample does not warn", {
  set.seed(11)
  y <- round(stats::rnorm(37), 1)
  psi <- function(theta) suppressWarnings(ee_percentile(theta, y = y, q = 0.5))
  expect_no_warning(
    estimate(MEstimator(stacked_equations = psi, init = stats::median(y)))
  )
})

test_that("the Rogan-Gladen misclassification fit does not warn", {
  d <- data.frame(
    r = c(1, 1, 0, 0, 0, 0),
    y = c(0, 0, 1, 1, 0, 0),
    y_star = c(1, 0, 1, 0, 1, 0),
    n = c(680, 270, 204, 38, 18, 71)
  )
  d <- d[rep(seq_len(nrow(d)), d$n), ]
  psi <- function(theta) {
    ee_rogan_gladen(theta, y = d$y, y_star = d$y_star, r = d$r)
  }
  expect_no_warning(
    estimate(MEstimator(stacked_equations = psi, init = rep(0.75, 4)))
  )
})

test_that("a stacked odds-ratio fit does not warn", {
  set.seed(42)
  n <- 500
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.8 * x))
  design <- cbind(1, x)
  psi <- function(theta) {
    rbind(
      ee_regression(theta[1:2], X = design, y = y, model = "logistic"),
      matrix(rep(exp(theta[2]) - theta[3], n), nrow = 1)
    )
  }
  expect_no_warning(
    estimate(MEstimator(stacked_equations = psi, init = c(0, 0, 1)))
  )
})

test_that("the regression calibration fit does not warn", {
  skip_if_not_installed("nleqslv")
  ref <- load_fixture("ee_regression_calibration")
  x_main <- cbind(1, ref$a_star)
  psi <- function(theta) {
    y_safe <- ifelse(ref$r == 0, -999, ref$y)
    ee_logit <- ee_regression(
      theta[4:5],
      X = x_main,
      y = y_safe,
      model = "logistic"
    )
    rbind(
      ee_regression_calibration(
        theta[1:3],
        beta = theta[5],
        a = ref$a,
        a_star = ref$a_star,
        r = ref$r
      ),
      ee_logit * rep(ref$r, each = nrow(ee_logit))
    )
  }
  expect_no_warning(
    estimate(
      MEstimator(stacked_equations = psi, init = ref$init),
      solver = "nleqslv"
    )
  )
})

# ---- nleqslv termination codes -----------------------------------------------
#
# nleqslv declares convergence through ftol, an absolute bound on the largest
# absolute function value with a default of 1e-8. A stack whose terms are large
# in magnitude cannot reach it however exactly the equations are solved, so such
# a fit stops on the x-value criterion instead and reports code 2, or stalls and
# reports code 3. The manual documents both as leaving the function values for
# the caller to judge rather than as failures, so deli judges them on a
# scale-free criterion. Weighting the logistic fixture by 1e12 puts the
# per-observation contributions at that magnitude while leaving the solution
# untouched.

weighted_logistic_psi <- function(weight) {
  ref <- load_fixture("ee_solver_lm_logistic")
  function(theta) {
    ee_regression(
      theta,
      X = ref$X,
      y = ref$y,
      model = "logistic",
      weights = rep(weight, length(ref$y))
    )
  }
}

test_that("nleqslv stays quiet when a solved fit cannot reach the absolute ftol", {
  skip_if_not_installed("nleqslv")
  ref <- load_fixture("ee_solver_lm_logistic")
  psi <- weighted_logistic_psi(1e12)
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = ref$theta),
      solver = "nleqslv"
    )
  })
  # The returned point is the solution, to fifteen digits.
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-12)
  # The absolute score is well above the floor an absolute check would use, so
  # only the scale-free criterion keeps this quiet.
  expect_gt(max(abs(rowSums(psi(m@theta)))), 1e-4)
})

test_that("nleqslv warns when a soft termination code leaves the equations unsolved", {
  skip_if_not_installed("nleqslv")
  # theta^2 + 1 has no real root, and nleqslv stalls rather than failing hard.
  psi <- function(theta) matrix(rep((theta[1]^2 + 1) / 10, 10), nrow = 1)
  expect_warning(
    estimate(
      MEstimator(stacked_equations = psi, init = 1),
      solver = "nleqslv"
    ),
    class = "deli_solver_not_converged"
  )
})

test_that("nleqslv warns on a hard termination code", {
  skip_if_not_installed("nleqslv")
  y <- rep(5, 40)
  psi <- function(theta) matrix(inverse_logit(theta[1]) - y, nrow = 1)
  expect_warning(
    estimate(
      MEstimator(stacked_equations = psi, init = 0),
      solver = "nleqslv"
    ),
    class = "deli_solver_not_converged"
  )
})

# rootSolve is the solver that returns a spurious root silently, so no
# non-convergence message may send a user to it as the remedy.
test_that("the nleqslv non-convergence message does not recommend rootSolve", {
  skip_if_not_installed("nleqslv")
  y <- rep(5, 40)
  psi <- function(theta) matrix(inverse_logit(theta[1]) - y, nrow = 1)
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = 0), solver = "nleqslv")
  )
  expect_length(seen, 1L)
  expect_no_match(conditionMessage(seen[[1]]), "rootSolve", fixed = TRUE)
})

test_that("the lm non-convergence message does not recommend rootSolve", {
  skip_if_not_installed("minpack.lm")
  ref <- load_fixture("ee_solver_lm_logistic")
  psi <- function(theta) {
    ee_regression(theta, X = ref$X, y = ref$y, model = "logistic")
  }
  seen <- collect_warnings(
    estimate(
      MEstimator(stacked_equations = psi, init = ref$init),
      solver = "lm",
      maxiter = 2
    )
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_no_match(conditionMessage(seen[[1]]), "rootSolve", fixed = TRUE)
})

# ---- minpack.lm termination codes --------------------------------------------
#
# Code 4 reports that the residual vector is orthogonal to every column of the
# Jacobian to within gtol, which is a stall rather than a solved system, so the
# returned point is judged rather than accepted.

test_that("the lm solver warns on a stall that leaves the equations unsolved", {
  skip_if_not_installed("minpack.lm")
  # theta^2 + 1 has no real root, and MINPACK stops on the orthogonality test.
  psi <- function(theta) matrix(rep((theta[1]^2 + 1) / 10, 10), nrow = 1)
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = 1), solver = "lm")
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(conditionMessage(seen[[1]]), "code 4", fixed = TRUE)
  expect_match(
    conditionMessage(seen[[1]]),
    "stopped without solving",
    fixed = TRUE
  )
})

test_that("the lm solver does not warn on a fit that meets a convergence test", {
  skip_if_not_installed("minpack.lm")
  ref <- load_fixture("ee_solver_lm_logistic")
  psi <- function(theta) {
    ee_regression(theta, X = ref$X, y = ref$y, model = "logistic")
  }
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = ref$init),
      solver = "lm"
    )
  })
  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-6)
})

# ---- GMM moment quality ------------------------------------------------------
#
# minimize_gmm only inspects optim's status code, and optim reports convergence
# on the flat tail of an objective that never reaches zero. A just-identified
# problem has as many moment conditions as parameters, so the moments must
# vanish at a solution, and the returned point is judged on the Newton step it
# would still take towards a root rather than on the size of its moments, which
# carries the scale of the estimating functions.
#
# The negative control below is the one that matters. Every stack built from a
# design matrix has mixed-sign per-observation contributions, so the
# non-cancellation share is small however wrong the parameters are, and any
# criterion built on that share cannot see the failure.

test_that("GMM warns when a just-identified fit does not solve its moments", {
  y <- rep(5, 40)
  psi <- function(theta) matrix(inverse_logit(theta[1]) - y, nrow = 1)
  expect_warning(
    estimate(GMMEstimator(stacked_equations = psi, init = 0)),
    class = "deli_solver_not_converged"
  )
})

# Five collinear mtcars columns. The default solver reports convergence, the
# intercept comes back at roughly a twentieth of its least-squares value, and
# every summed moment stays small because the residual contributions cancel.
collinear_mtcars_psi <- function() {
  design <- cbind(1, mtcars$wt, mtcars$hp, mtcars$disp, mtcars$drat)
  function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
}

test_that("GMM warns when a design-matrix fit stops far from the solution", {
  psi <- collinear_mtcars_psi()
  seen <- collect_warnings({
    m <- estimate(GMMEstimator(stacked_equations = psi, init = rep(0, 5)))
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(conditionMessage(seen[[1]]), "Newton step", fixed = TRUE)
  # The returned intercept is nowhere near the least-squares one.
  ols <- stats::coef(stats::lm(mpg ~ wt + hp + disp + drat, data = mtcars))
  expect_gt(abs(m@theta[[1]] - ols[[1]]), 20)
})

test_that("is_root() cannot see a design-matrix fit that stopped short", {
  psi <- collinear_mtcars_psi()
  suppressWarnings({
    m <- estimate(GMMEstimator(stacked_equations = psi, init = rep(0, 5)))
  })
  evald <- psi(m@theta)
  # is_root() is the criterion the M-estimation solvers are judged on. It passes
  # here, which is why the moments need a scale-aware criterion of their own.
  expect_true(is_root(evald, m@theta))
  expect_lt(
    max(abs(rowSums(evald)) / rowSums(abs(evald))),
    noncancel_ceiling
  )
  # Nor can the absolute size of the moments be the signal. The largest summed
  # moment here is under 5, while the weighted logistic fit above is solved to
  # fifteen digits with a summed score of 1e7, so no absolute threshold
  # separates the two.
  expect_lt(max(abs(rowSums(evald))), 5)
})

test_that("a just-identified GMM fit that solves its moments does not warn", {
  design <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
  expect_no_warning({
    m <- estimate(GMMEstimator(stacked_equations = psi, init = c(0, 0, 0)))
  })
  expect_equal(
    unname(m@theta),
    unname(stats::coef(stats::lm(mpg ~ wt + hp, data = mtcars))),
    tolerance = 1e-6
  )
})

# An over-identified problem cannot drive every moment to zero, and a subset fit
# holds parameters at their initial values while still summing every equation
# into the objective. Neither can be judged by the size of its moments.
test_that("an over-identified GMM fit is not judged by the size of its moments", {
  set.seed(5)
  y <- stats::rnorm(200, mean = 3, sd = sqrt(2))
  # Two moment conditions for one parameter: the variance condition cannot hold
  # exactly at the value the mean condition picks out.
  psi <- function(theta) {
    rbind(y - theta[1], (y - theta[1])^2 - 2)
  }
  expect_no_warning({
    m <- estimate(GMMEstimator(stacked_equations = psi, init = 2))
  })
  expect_gt(max(abs(rowSums(psi(m@theta)))), 1)
})

test_that("a subset M fit is judged only on the equations it solves", {
  set.seed(19)
  y <- stats::rnorm(200, mean = 5)
  psi <- function(theta) {
    rbind(y - theta[1], (y - theta[1])^2 - theta[2])
  }
  # theta[1] is held at zero, far from the mean of y, so the first equation
  # cannot vanish however well the fit went. Only the second is solved, so only
  # the second may be judged.
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 1), subset = 2L)
    )
  })
  expect_equal(unname(m@theta), c(0, mean(y^2)), tolerance = 1e-6)
  # Judging the whole stack instead would report the first equation.
  theta <- unname(m@theta)
  bread <- compute_bread(psi, theta, "capprox", 1e-9) / 200
  expect_equal(unsolved_point(psi(theta), theta, bread)$row, 1L)
})

test_that("a subset GMM fit is not judged by the size of its moments", {
  design <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
  expect_no_warning({
    m <- estimate(
      GMMEstimator(
        stacked_equations = psi,
        init = c(0, 0, 0),
        subset = c(1L, 2L)
      )
    )
  })
  # The third parameter never moves, so the third moment cannot vanish.
  expect_gt(max(abs(rowSums(psi(m@theta)))), 1)
})

# A non-differentiable moment condition has a singular bread everywhere, so the
# distance to a root cannot be measured with it and the point stays unjudged.
# This is the excuse the M-estimation path already makes for the same equations.
test_that("a non-differentiable GMM fit is left unjudged", {
  set.seed(11)
  y <- round(stats::rnorm(37), 1)
  psi <- function(theta) suppressWarnings(ee_percentile(theta, y = y, q = 0.5))
  expect_no_warning({
    m <- suppressWarnings(estimate(
      GMMEstimator(stacked_equations = psi, init = stats::median(y))
    ))
  })
  # Every contribution to the median condition is 0.5 in magnitude, so an
  # odd-sized sample cannot sum to zero however well the fit went.
  expect_gte(max(abs(rowSums(psi(m@theta)))), 0.5)
})

test_that("a GMM fit whose minimizer failed is not warned about twice", {
  design <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
  seen <- collect_warnings(
    estimate(
      GMMEstimator(stacked_equations = psi, init = c(0, 0, 0)),
      maxiter = 2
    )
  )
  # optim reports its own failure, so the moments are not judged on top of it.
  expect_length(seen, 1L)
  expect_match(conditionMessage(seen[[1]]), "did not converge")
})

test_that("a custom GMM minimizer has its moments judged", {
  y <- rep(5, 40)
  psi <- function(theta) matrix(inverse_logit(theta[1]) - y, nrow = 1)
  # A minimizer that returns the starting values and reports nothing.
  stalled <- function(stacked_equations, init) init
  expect_warning(
    estimate(
      GMMEstimator(stacked_equations = psi, init = 0),
      solver = stalled
    ),
    class = "deli_solver_not_converged"
  )
})

# ---- equations with no variation across observations --------------------------
#
# An equation whose contributions are the same at every observation is a
# documented way of writing a relation among the parameters: the ratio of two
# means in vignette("custom-estimating-equations") and the causal effect rows of
# ee_ipw, ee_aipw and ee_gformula are all built that way. The non-cancellation
# share of such a row is exactly one at every point except an exact zero, so it
# says nothing about how close the point is to a root and cannot be the criterion
# for the row. What must vanish there is the repeated value itself, and that is
# what unsolved_equation() measures, against the larger of the magnitude of the
# estimates and the largest contribution anywhere in the stack. Taking the larger
# of the two is what frees the reading from the scale of the estimating
# functions.

# The ratio of two means, exactly as vignette("custom-estimating-equations")
# writes it: an equation for each mean and a constant third row for the ratio.
# `scale` multiplies the whole stack, which changes nothing about the solution.
ratio_psi <- function(seed, n = 200, scale = 1) {
  set.seed(seed)
  y1 <- stats::rnorm(n, mean = 4, sd = 1)
  y2 <- stats::rnorm(n, mean = 2, sd = 1)
  function(theta) {
    rbind(
      scale * (y1 - theta[1]),
      scale * (y2 - theta[2]),
      rep(scale * (theta[1] / theta[2] - theta[3]), n)
    )
  }
}

test_that("a GMM fit of the ratio stack does not warn at any seed", {
  # Correct to seven digits against the M-estimation solution at every seed, so
  # any warning here is a false alarm rather than a fit worth reporting.
  for (seed in 1:40) {
    psi <- ratio_psi(seed)
    expect_no_warning({
      g <- estimate(GMMEstimator(stacked_equations = psi, init = c(1, 1, 1)))
    })
    expect_no_warning({
      m <- estimate(MEstimator(stacked_equations = psi, init = c(1, 1, 1)))
    })
    expect_equal(unname(g@theta), unname(m@theta), tolerance = 1e-5)
  }
})

test_that("the causal estimators do not warn when refitted as GMM", {
  set.seed(42)
  n <- 400
  w1 <- stats::rnorm(n)
  w2 <- stats::rbinom(n, 1, 0.4)
  a <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.5 * w1 + 0.3 * w2))
  y <- 2 + 1.5 * a + w1 - 0.5 * w2 + stats::rnorm(n)
  w_ps <- cbind(1, w1, w2)
  x <- cbind(1, a, w1, w2)
  x1 <- cbind(1, 1, w1, w2)
  x0 <- cbind(1, 0, w1, w2)
  cases <- list(
    list(function(theta) ee_ipw(theta, y = y, A = a, W = w_ps), rep(0, 6)),
    list(
      function(theta) ee_gformula(theta, y = y, X = x, X1 = x1, X0 = x0),
      rep(0, 7)
    ),
    list(
      function(theta) {
        ee_aipw(theta, y = y, A = a, W = w_ps, X = x, X1 = x1, X0 = x0)
      },
      rep(0, 10)
    )
  )
  for (case in cases) {
    expect_no_warning(
      estimate(GMMEstimator(stacked_equations = case[[1]], init = case[[2]]))
    )
  }
})

test_that("is_root() reaches the same verdict when the equations are rescaled", {
  y <- c(3, 5, 7, 9)
  theta <- c(6, 6 - 1e-9)
  # A mean equation solved exactly and a constant contrast solved to within a
  # billionth, which is a root by any reasonable reading.
  stack <- function(scale) {
    rbind(scale * (y - theta[1]), rep(scale * (theta[1] - theta[2]), 4))
  }
  for (scale in c(1e6, 1e9, 1e12)) {
    # Above the absolute floor at every scale used, so the floor is not what
    # decides these.
    expect_gt(max(abs(rowSums(stack(scale)))), score_floor)
    expect_true(is_root(stack(scale), theta))
  }
})

test_that("is_root() rejects an equation with no variation that does not vanish", {
  y <- c(3, 5, 7, 9)
  theta <- c(6, 5)
  # The contrast row demands theta[1] == theta[2] and is out by a whole unit.
  stack <- function(scale) {
    rbind(scale * (y - theta[1]), rep(scale * (theta[1] - theta[2]), 4))
  }
  for (scale in c(1, 1e6, 1e12)) {
    expect_false(is_root(stack(scale), theta))
  }
})

test_that("unsolved_equation() names the equation that failed", {
  # Row 1's contributions all point the same way, so it cannot cancel and is
  # judged on its size against the scale of the stack, which it misses. Row 2
  # leaves a summed score three hundred times larger but cancels well, so it
  # passes. The number reported must be row 1's, not the largest score in the
  # stack.
  ef <- rbind(
    c(100, 200, 150, 300, 250),
    c(1e6, -1e6, 1e6, -1e6, 3e5)
  )
  found <- unsolved_equation(ef, c(1, 1))
  expect_equal(found$row, 1L)
  expect_false(found$constant)
  expect_equal(found$score, 1000)
  expect_gt(max(abs(rowSums(ef))), 1e5)
})

test_that("the lm solver stays quiet on a rescaled stack it solves exactly", {
  skip_if_not_installed("minpack.lm")
  # Multiplying the stack by 1e10 puts the constant row's round-off residual
  # above the absolute score floor, which is where the cancellation share alone
  # used to report a failure on an exactly solved fit. minpack reports a soft
  # code here only some of the time, so the point is judged directly as well.
  psi <- ratio_psi(13, scale = 1e10)
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(1, 1, 1)),
      solver = "lm"
    )
  })
  set.seed(13)
  y1 <- stats::rnorm(200, mean = 4, sd = 1)
  y2 <- stats::rnorm(200, mean = 2, sd = 1)
  expect_equal(
    unname(m@theta),
    c(mean(y1), mean(y2), mean(y1) / mean(y2)),
    tolerance = 1e-12
  )
  evald <- psi(m@theta)
  expect_gt(max(abs(rowSums(evald))), score_floor)
  expect_true(is_root(evald, m@theta))
})

test_that("GMM warns when a constant moment condition cannot vanish", {
  # The third equation demands that a positive ratio equal a negative number, so
  # no parameter value solves it.
  set.seed(11)
  a <- stats::rnorm(200, mean = 4)
  b <- stats::rnorm(200, mean = 2)
  psi <- function(theta) {
    rbind(
      a - theta[1],
      b - theta[2],
      rep(theta[1] / theta[2] + theta[3]^2 + 1, 200)
    )
  }
  seen <- collect_warnings(
    estimate(GMMEstimator(stacked_equations = psi, init = c(1, 1, 1)))
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(
    conditionMessage(seen[[1]]),
    "Moment condition 3 has the same value at every observation",
    fixed = TRUE
  )
})

# The share is the only thing that sees a gross failure on a stack with no
# constant row, so splitting the criterion must not cost it. A tobit fit refitted
# as GMM returns a log scale parameter of 10.6 against -0.002 and a slope of 1.37
# against 0.478, at a Newton step well under the ceiling.
test_that("GMM still warns on a tobit fit that stops far from the solution", {
  set.seed(123)
  x <- cbind(1, stats::rnorm(200))
  y <- pmax(1 + 0.5 * x[, 2] + stats::rnorm(200), 0)
  psi <- function(theta) ee_tobit(theta, X = x, y = y, lower = 0)
  init <- c(mean(y), 0, log(stats::sd(y)))
  seen <- collect_warnings({
    g <- estimate(GMMEstimator(stacked_equations = psi, init = init))
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(conditionMessage(seen[[1]]), "do not cancel", fixed = TRUE)
  m <- estimate(MEstimator(stacked_equations = psi, init = init))
  expect_gt(abs(g@theta[[3]] - m@theta[[3]]), 10)
  # The Newton step cannot see this one, which is why the share is kept.
  bread <- compute_bread(psi, g@theta, "capprox", 1e-9) / 200
  step <- relative_newton_step(bread, rowSums(psi(g@theta)) / 200, g@theta)
  expect_lt(step, newton_step_ceiling)
})

# ---- equations whose contributions all carry one sign -------------------------
#
# The share of an equation's contributions that fails to cancel is exactly one
# whenever they all carry the same sign, whatever else varies across
# observations, so for such an equation the share says nothing about how close
# the point is to a root. An equation with no variation is only one way to reach
# that state. The relation rows of the ratio stack and of the causal estimators
# multiplied by an observation weight are another: they are not n copies of one
# value, so an exact-equality test does not recognise them, yet they cannot
# cancel and their share is one at a point that solves them to the last
# representable digit.

# The ratio of two weighted means, written as the ratio stack of
# vignette("custom-estimating-equations") with an observation weight on every
# row, including the relation row.
weighted_ratio_psi <- function(seed, n = 200, scale = 1e10) {
  set.seed(seed)
  y1 <- stats::rnorm(n, mean = 4, sd = 1)
  y2 <- stats::rnorm(n, mean = 2, sd = 1)
  wt <- stats::runif(n, 0.5, 2)
  function(theta) {
    scale *
      rbind(
        wt * (y1 - theta[1]),
        wt * (y2 - theta[2]),
        wt * (theta[1] / theta[2] - theta[3])
      )
  }
}

# The solution of that stack in closed form.
weighted_ratio_theta <- function(seed, n = 200) {
  set.seed(seed)
  y1 <- stats::rnorm(n, mean = 4, sd = 1)
  y2 <- stats::rnorm(n, mean = 2, sd = 1)
  wt <- stats::runif(n, 0.5, 2)
  means <- c(sum(wt * y1), sum(wt * y2)) / sum(wt)
  c(means, means[[1]] / means[[2]])
}

test_that("is_root() judges a one-sided equation by its size, not by its share", {
  y <- c(3, 5, 7, 9)
  wt <- c(0.7, 1.3, 0.9, 1.1)
  theta <- c(6, 6 - 1e-9)
  # A mean equation solved exactly and a weighted contrast solved to within a
  # billionth, which is a root by any reasonable reading. The contrast row takes
  # a different value at every observation, so no exact-equality test sees it,
  # and every contribution carries the same sign, so its share is exactly one.
  stack <- function(scale) {
    rbind(scale * (y - theta[1]), scale * wt * (theta[1] - theta[2]))
  }
  for (scale in c(1e6, 1e9, 1e12)) {
    ef <- stack(scale)
    expect_gt(max(abs(rowSums(ef))), score_floor)
    expect_false(isTRUE(all(ef[2, ] == ef[2, 1L])))
    expect_equal(abs(sum(ef[2, ])) / sum(abs(ef[2, ])), 1)
    expect_true(is_root(ef, theta))
  }
})

test_that("is_root() rejects a one-sided equation that does not vanish", {
  y <- c(3, 5, 7, 9)
  wt <- c(0.7, 1.3, 0.9, 1.1)
  theta <- c(6, 5)
  # The weighted contrast demands theta[1] == theta[2] and is out by a unit.
  stack <- function(scale) {
    rbind(scale * (y - theta[1]), scale * wt * (theta[1] - theta[2]))
  }
  for (scale in c(1, 1e6, 1e12)) {
    expect_false(is_root(stack(scale), theta))
  }
})

test_that("a weighted ratio stack solved to the last digit does not warn", {
  skip_if_not_installed("nleqslv")
  for (seed in 1:6) {
    for (solver in c("rootSolve", "nleqslv")) {
      expect_no_warning({
        m <- estimate(
          MEstimator(
            stacked_equations = weighted_ratio_psi(seed),
            init = c(1, 1, 1)
          ),
          solver = solver
        )
      })
      expect_equal(
        unname(m@theta),
        weighted_ratio_theta(seed),
        tolerance = 1e-12
      )
    }
  }
})

test_that("a weighted ratio stack is judged when its relation does not vanish", {
  # The same stack with a relation no parameter value can satisfy: a ratio of
  # two positive means cannot equal minus one times a square plus one.
  set.seed(4)
  n <- 200
  y1 <- stats::rnorm(n, mean = 4, sd = 1)
  y2 <- stats::rnorm(n, mean = 2, sd = 1)
  wt <- stats::runif(n, 0.5, 2)
  psi <- function(theta) {
    rbind(
      wt * (y1 - theta[1]),
      wt * (y2 - theta[2]),
      wt * (theta[1] / theta[2] + theta[3]^2 + 1)
    )
  }
  seen <- collect_warnings(
    estimate(GMMEstimator(stacked_equations = psi, init = c(1, 1, 1)))
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  # The relation row is one-signed but not constant, so it is judged on its mean
  # against the scale of the problem. The message must name that reading rather
  # than the share reading, which did not judge this row.
  expect_match(
    conditionMessage(seen[[1]]),
    "contributions of moment condition 3 all carry one sign",
    fixed = TRUE
  )
  expect_no_match(conditionMessage(seen[[1]]), "do not cancel", fixed = TRUE)
})

# ---- solver reports that are judged rather than trusted -----------------------
#
# Every convergence test a solver applies is relative to something of its own:
# the residuals it started from, the step it last took, or an absolute bound on
# the function values. None of them is a statement that the estimating equations
# are solved, so a report of success is judged against the returned point rather
# than accepted. The judgement is made after the bread, because the Jacobian is
# what measures the distance still to travel and it is computed for the sandwich
# regardless.

test_that("the lm solver warns when a met convergence test leaves the equations unsolved", {
  skip_if_not_installed("minpack.lm")
  set.seed(42)
  n <- 400
  w1 <- stats::rnorm(n)
  w2 <- stats::rbinom(n, 1, 0.4)
  a <- stats::rbinom(n, 1, stats::plogis(-0.5 + 0.5 * w1 + 0.3 * w2))
  y <- 2 + 1.5 * a + w1 - 0.5 * w2 + stats::rnorm(n)
  x <- cbind(1, a, w1, w2)
  x1 <- cbind(1, 1, w1, w2)
  x0 <- cbind(1, 0, w1, w2)
  # Multiplying by 1e8 leaves the solution untouched and puts the residuals
  # where nls.lm reports its relative sum-of-squares test met at the starting
  # values, which is code 1 rather than the stall code 4.
  psi <- function(theta) {
    1e8 * ee_gformula(theta, y = y, X = x, X1 = x1, X0 = x0)
  }
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = rep(0, 7)),
      solver = "lm"
    )
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  # It never left the starting values, against a solution near (1.5, 1.9, 1.5).
  expect_lt(max(abs(m@theta)), 1e-6)
})

test_that("nleqslv is judged when its function criterion is met without moving", {
  skip_if_not_installed("nleqslv")
  set.seed(4)
  y <- stats::rnorm(200, mean = 5)
  # Scaling a mean equation down by 1e-12 puts every function value under
  # nleqslv's ftol, an absolute bound with a default of 1e-8, at the starting
  # values. nleqslv therefore reports code 1, its function criterion met,
  # without moving at all. The code is a statement about the size of the
  # function values rather than about the parameters, so the point is judged.
  psi <- function(theta) matrix(1e-12 * (y - theta[1]), nrow = 1)
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = 0),
      solver = "nleqslv"
    )
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(conditionMessage(seen[[1]]), "(code 1)", fixed = TRUE)
  # It never left zero, against a solution of 5.01.
  expect_equal(unname(m@theta), 0)
  expect_gt(mean(y), 5)
  # The contributions alone cannot see this one: the summed score is 1e-9, well
  # under the floor. Only the Newton step measures the distance still to travel.
  expect_true(is_root(psi(0), 0))
})

test_that("a custom M solver has its returned point judged", {
  set.seed(7)
  w <- stats::rnorm(200)
  y <- 1 + 2 * w + stats::rnorm(200)
  # A stack whose third row demands theta[3] == 7, handed to a solver that
  # returns the starting values and reports nothing at all.
  psi <- function(theta) {
    rbind(
      ee_regression(theta[1:2], X = cbind(1, w), y = y, model = "linear"),
      rep(theta[3] - 7, length(y))
    )
  }
  stalled <- function(stacked_equations, init) init
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0, 0)),
      solver = stalled
    )
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(
    conditionMessage(seen[[1]]),
    "The solver returned values that do not solve the estimating equations",
    fixed = TRUE
  )
  # No built-in solver is named, because none was used, and no message may send
  # a user to rootSolve in any case.
  expect_no_match(conditionMessage(seen[[1]]), "rootSolve", fixed = TRUE)
  expect_equal(unname(m@theta), c(0, 0, 0))
})

test_that("a custom M solver that solves the equations stays quiet", {
  set.seed(7)
  w <- stats::rnorm(200)
  y <- 1 + 2 * w + stats::rnorm(200)
  psi <- function(theta) {
    rbind(
      ee_regression(theta[1:2], X = cbind(1, w), y = y, model = "linear"),
      rep(theta[3] - 7, length(y))
    )
  }
  # The same stack and a solver that also reports nothing, but reaches the
  # root. Judging a custom solver's point must not cost this fit its silence.
  minimizing <- function(stacked_equations, init) {
    stats::optim(
      init,
      function(theta) sum(stacked_equations(theta)^2),
      method = "BFGS",
      control = list(reltol = 1e-14, maxit = 2000)
    )$par
  }
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0, 0)),
      solver = minimizing
    )
  })
  expect_equal(unname(m@theta)[[3]], 7)
})

test_that("rootSolve warns when it never leaves starting values that do not solve the equations", {
  set.seed(42)
  n <- 200
  w <- stats::rnorm(n)
  # A response large enough to defeat multiroot's own Jacobian, centred so that
  # every row's contributions are mixed-sign and the non-cancellation share
  # cannot see the failure. multiroot returns the starting values, which is the
  # one case its own report of a failed convergence test cannot be taken at face
  # value, because a non-differentiable equation legitimately returns them.
  #
  # deriv_method = "exact" is load-bearing rather than incidental. At this data
  # scale the default finite-difference bread collapses to a matrix of zeros,
  # dx being far smaller than the resolution of the contributions, and a bread
  # that cannot be solved leaves the point unjudged. The default path therefore
  # does not catch this fit.
  y <- 1e9 * (0.7 * w + stats::rnorm(n))
  x <- cbind(1, w)
  psi <- function(theta) ee_regression(theta, X = x, y = y, model = "linear")
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0)),
      deriv_method = "exact"
    )
  })
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_equal(unname(m@theta), c(0, 0))
  # Neither reading of the contributions alone can see this one.
  expect_true(is_root(psi(m@theta), unname(m@theta)))
})

test_that("a fit whose Jacobian does not exist is left unjudged", {
  set.seed(2024)
  y <- 0.6 + stats::rnorm(200)
  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }
  init <- c(1, stats::median(y))
  # multiroot reports its convergence test failing here and returns the starting
  # values, which is the right answer for a median. The bread is singular, so the
  # distance to a root cannot be measured with it and the point stays unjudged.
  expect_no_warning({
    m <- estimate(MEstimator(stacked_equations = psi, init = init))
  })
  expect_equal(unname(m@theta), unname(init))
  expect_true(is.na(relative_newton_step(
    m@bread,
    rowSums(psi(unname(m@theta))) / 200,
    unname(m@theta)
  )))
})

test_that("a moment condition that does not sum to a finite value is named as such", {
  # A minimizer that stops where the first moment condition is undefined. It
  # reports nothing, so its point is judged like one a built-in minimizer
  # declared solved.
  psi <- function(theta) {
    matrix(c(1 / (theta[1] - 2), rep(theta[1] - 1, 9)), nrow = 1)
  }
  stalled <- function(stacked_equations, init) 2
  seen <- collect_warnings(
    estimate(GMMEstimator(stacked_equations = psi, init = 0), solver = stalled)
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(conditionMessage(seen[[1]]), "is not finite", fixed = TRUE)
  expect_no_match(conditionMessage(seen[[1]]), "do not cancel", fixed = TRUE)
})

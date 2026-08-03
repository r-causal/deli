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
  reported <- reported_call(err)

  expect_match(reported, "method(estimate, deli::MEstimator)", fixed = TRUE)
  expect_false(grepl("estimate_m_estimator", reported, fixed = TRUE))
})

# The GMM path has the same shape and needs the same treatment: the refusal of a
# solver that is neither a name nor a function is raised inside the minimization
# dispatcher, two frames below the method the caller reached.
test_that("an invalid GMM solver type reports a frame the user can name", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    matrix(y - theta[1], nrow = 1)
  }

  g <- GMMEstimator(stacked_equations = psi, init = c(0))
  err <- tryCatch(estimate(g, solver = 42), error = function(e) e)
  reported <- reported_call(err)

  expect_match(conditionMessage(err), "must be a string or function")
  expect_match(reported, "method(estimate, deli::GMMEstimator)", fixed = TRUE)
  expect_false(grepl("minimize_gmm", reported, fixed = TRUE))
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
# These tests pin the API surface and numerical behavior of the
# Levenberg-Marquardt solver. In Python delicatessen it is the default,
# reached through scipy.optimize.root(method = "lm"); the deli surface
# exposes it as solver = "lm".
#
# Fixtures ee_solver_lm_linear and ee_solver_lm_logistic come from Python
# solver = "lm" (fixtures/generate_solver_lm_fixtures.py). The linear
# problem is well-conditioned so every root-finder agrees to tight
# tolerance; the logistic problem is a nonlinear system that exercises the
# least-squares behavior of Levenberg-Marquardt.

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
# itself emits a "steady-state not reached" warning, which is muffled where it is
# raised and reworded in deli's own terms, and an exhausted budget is not what
# that report says, so it has to be surfaced from the iteration count separately.
# The lm and nleqslv branches already warn on their own non-convergence codes.
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
# warning rootSolve raises alongside it is muffled where it is raised, so what
# the solve reports rests on that warning being read rather than only
# suppressed. A beta regression started from a scale parameter that is far too
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
  # The median equation cannot be searched, so the fit holds wherever it is
  # started and both parameters have to be started at the values the data give
  # them. The deviation equation is solved at twice the mean deviation above the
  # median, and its residual there is round-off while the median equation's is
  # not, which is the pair this block exists to keep quiet. Starting the
  # deviation anywhere else returns a value nobody solved for; see "a singular
  # bread does not excuse the equation whose own row is not flat" below.
  init <- c(2 * mean(pmax(y - stats::median(y), 0)), stats::median(y))
  expect_no_warning(
    estimate(MEstimator(stacked_equations = psi, init = init))
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

# The share is what names the equation on a stack with no constant row, so
# splitting the criterion must not cost it. A tobit fit refitted as GMM returns a
# log scale parameter of 10.6 against -0.002 and a slope of 1.37 against 0.478,
# and the report has to say which moment condition that is.
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
  # The share is asked before the Newton step and names a row, so the step is
  # never taken here however far past the ceiling it is. That ordering is what
  # keeps the report naming the equation rather than quoting a step.
  bread <- compute_bread(psi, g@theta, "capprox", 1e-9) / 200
  step <- relative_newton_step(bread, rowSums(psi(g@theta)) / 200, g@theta)
  expect_gt(step, newton_step_ceiling)
  expect_false(is.na(unsolved_point(psi(g@theta), g@theta, bread)$row))
})

# ---- equations whose contributions all carry one sign -------------------------
#
# The share of an equation's contributions that fails to cancel is exactly one
# whenever they all carry the same sign, whatever else varies across
# observations, so for such an equation the share says nothing about how close
# the point is to a root. An equation with no variation is only one way to reach
# that state. The relation rows of the ratio stack and of the causal estimators
# multiplied by an observation weight are another: they are not n copies of one
# value, so an exact-equality test does not recognize them, yet they cannot
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
# than accepted. The judgment is made after the bread, because the Jacobian is
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
  # A response large enough to defeat multiroot's own Jacobian, centered so that
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
  init <- c(2 * mean(pmax(y - stats::median(y), 0)), stats::median(y))
  # multiroot reports its convergence test failing here and returns the starting
  # values, which is the right answer for a median. The bread is singular, so no
  # Newton step can be solved for and that reading says nothing, which is what is
  # pinned below. The reading that does judge the point is the step
  # relative_singular_step() takes along the directions the bread still sees, and
  # it accepts it: both equations are solved where the fit stopped.
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

test_that("a singular bread does not excuse the equation whose own row is not flat", {
  set.seed(2024)
  y <- 0.6 + stats::rnorm(200)
  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }
  # The median equation cannot be searched, so multiroot returns the starting
  # values with its convergence test failing. That leaves the first parameter
  # where the caller put it, and the caller put it a long way from the value the
  # data give it: the deviation equation reads 2 * mean((y - median)+), which is
  # 0.79 against the 1 it was started at.
  init <- c(1, stats::median(y))
  seen <- collect_warnings({
    m <- estimate(MEstimator(stacked_equations = psi, init = init))
  })
  theta <- unname(coef(m))
  expect_equal(theta, unname(init))
  expect_equal(2 * mean(pmax(y - theta[[2]], 0)), 0.790488, tolerance = 1e-5)
  # No other reading sees it. The contributions cancel well enough to say
  # nothing, the flat row belongs to the median equation, which is solved where
  # the fit stopped, and the whole bread cannot be solved for a Newton step.
  ef <- psi(theta)
  expect_true(is_root(ef, theta))
  expect_equal(rowSums(ef)[[2]], 0)
  expect_true(all(m@bread[2, ] == 0))
  expect_true(is.na(relative_newton_step(
    unname(m@bread),
    rowSums(ef) / 200,
    theta
  )))
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
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

# ---- a bread lost to the scale of the estimating functions --------------------
#
# The step floor in approx_differentiation() rescues a step lost against a large
# `theta`. It cannot reach the fit below, where the step is applied exactly and
# the summed equations are the large quantity: every difference falls below their
# floating-point resolution, so the bread is a matrix of zeros and the variance
# is zero with it. Nothing about the step is wrong, so the loss has to be
# reported from the significance of the differences themselves. See the section
# of test-derivative.R that pins the same reading on the derivative alone.

test_that("a bread lost to the scale of the estimating functions is reported", {
  set.seed(42)
  n <- 200
  w <- stats::rnorm(n)
  # The same stack as the divergence test above, which reaches its own reading
  # only under `deriv_method = "exact"`. Under the default finite difference the
  # summed equations run to 1e11, whose neighboring doubles are 1.53e-05 away,
  # and the default step moves them by 4e-07.
  y <- 1e9 * (0.7 * w + stats::rnorm(n))
  x <- cbind(1, w)
  psi <- function(theta) ee_regression(theta, X = x, y = y, model = "linear")
  seen <- collect_warnings({
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0)))
  })
  lost <- Filter(
    function(cnd) inherits(cnd, "deli_finite_difference_lost"),
    seen
  )
  # One report for the operation, whatever the number of entries lost, which is
  # what without_repeated_warnings() exists to deliver.
  expect_length(lost, 1L)
  # The derivative is nowhere near zero, so neither is the loss a property of the
  # problem: exact differentiation at the same point returns a bread whose
  # diagonal is of order one.
  exact_bread <- compute_bread(psi, unname(m@theta), "exact") / n
  expect_gt(min(abs(diag(exact_bread))), 0.9)
})

test_that("a well-scaled fit reports no lost bread", {
  set.seed(42)
  n <- 200
  w <- stats::rnorm(n)
  # The same stack on data of an ordinary magnitude. The summed equations vanish
  # at the solution, so the differences taken there carry every digit they have
  # and nothing may be reported.
  y <- 0.7 * w + stats::rnorm(n)
  x <- cbind(1, w)
  psi <- function(theta) ee_regression(theta, X = x, y = y, model = "linear")
  expect_no_warning({
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0)))
  })
  expect_equal(
    unname(coef(m)),
    unname(stats::coef(stats::lm(y ~ w))),
    tolerance = 1e-6
  )
})

# ---- points a zero bread row leaves unjudged ---------------------------------
#
# relative_newton_step() is the only reading that sees a stack whose mixed-sign
# contributions cancel well at points that are not roots, and it needs a bread
# that can be solved. A bread carrying an identically zero row cannot be, so the
# exact step is `NA` and what is left is the contribution readings and the step
# relative_singular_step() takes along the directions the bread does see. Where
# the equation that owns the zero row is the one left unsolved, neither of those
# reaches it, and a wrong root would be returned in silence.
#
# The excuse a singular bread earns is owed to the equation that is solved at the
# point, not to the one that is not: the fit in "a fit whose Jacobian does not
# exist is left unjudged" above returns both equations' own roots and stays
# quiet, and must go on doing so.

test_that("a zero bread row does not excuse an equation that is not solved", {
  skip_if_not_installed("minpack.lm")
  set.seed(1)
  y <- stats::rnorm(50)
  seen <- collect_warnings({
    m <- m_estimate(
      function(t) ee_positive_mean_deviation(t, y),
      init = c(0, 0),
      solver = "lm"
    )
  })
  theta <- unname(coef(m))
  ef <- suppressWarnings(ee_positive_mean_deviation(theta, y))
  # The median equation contributes 0.5 - (y <= theta[2]) per observation, so
  # its summed score counts the observations either side of the returned value:
  # 37 of the 50 fall below it, leaving 25 - 37. The returned value is not the
  # sample median, 0.129, and sits between the 37th and 38th order statistics,
  # 0.697 and 0.738.
  expect_equal(rowSums(ef)[[2]], -12)
  expect_gt(abs(theta[[2]] - stats::median(y)), 0.5)
  # Neither contribution reading sees it: the contributions are mixed-sign and
  # 48% of their mass fails to cancel, under the 90% ceiling.
  expect_true(is_root(ef, theta))
  # And the Newton step cannot be taken, because the median equation's row of the
  # bread is identically zero.
  expect_true(all(m@bread[2, ] == 0))
  expect_true(is.na(relative_newton_step(
    unname(m@bread),
    rowSums(ef) / 50,
    theta
  )))
  # None of which may leave the point passing for a root.
  not_converged <- Filter(
    function(cnd) inherits(cnd, "deli_solver_not_converged"),
    seen
  )
  expect_length(not_converged, 1L)
  # The solver reports a convergence test that was met, so nothing it says names
  # the equation it left behind. The report has to name it and the score it left
  # there, which is the whole of what this reading contributes over the solver's
  # own account of itself.
  reported <- conditionMessage(not_converged[[1]])
  expect_match(
    reported,
    "Estimating equation 2 does not move when any parameter does",
    fixed = TRUE
  )
  expect_match(reported, "it sums to -12 rather than to zero", fixed = TRUE)
})

test_that("the non-differentiability warning still surfaces once beside the report", {
  skip_if_not_installed("minpack.lm")
  set.seed(1)
  y <- stats::rnorm(50)
  # One fit evaluates the estimating function many times and every evaluation
  # warns that the median is not differentiable, so the scope in R/conditions.R
  # has to go on delivering that warning once for the operation while the report
  # of the unsolved point arrives beside it.
  seen <- collect_warnings(
    m_estimate(
      function(t) ee_positive_mean_deviation(t, y),
      init = c(0, 0),
      solver = "lm"
    )
  )
  messages <- vapply(seen, conditionMessage, character(1))
  expect_length(grep("not differentiable", messages, fixed = TRUE), 1L)
  expect_length(
    Filter(function(cnd) inherits(cnd, "deli_solver_not_converged"), seen),
    1L
  )
  # The two of them and nothing else. A bread of zeros is not a bread lost to
  # rounding, so the reading of the section above may not fire here.
  expect_length(seen, 2L)
})

test_that("a fit whose bread has no zero row stays quiet under the lm solver", {
  skip_if_not_installed("minpack.lm")
  set.seed(11)
  y <- stats::rnorm(100, mean = 2)
  psi <- function(theta) ee_mean(theta, y = y)
  expect_no_warning({
    m <- m_estimate(stacked_equations = psi, init = 0, solver = "lm")
  })
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-8)
})

test_that("GMM names the moment condition its bread has gone flat under", {
  set.seed(1)
  y <- stats::rnorm(50)
  # The same stack as the M-estimation fits above, minimized rather than solved.
  # A flat row reaches the moment conditions the same way it reaches the
  # estimating equations, so the reading is shared and each path words its own
  # report.
  psi <- function(theta) ee_positive_mean_deviation(theta, y = y)
  seen <- collect_warnings({
    g <- estimate(GMMEstimator(stacked_equations = psi, init = c(0, 0)))
  })
  theta <- unname(coef(g))
  ef <- suppressWarnings(ee_positive_mean_deviation(theta, y))
  # 28 of the 50 observations fall below the returned value, leaving 25 - 28,
  # and the median moment's row of the bread is identically zero, so the
  # minimizer had no reading of it at all.
  expect_equal(rowSums(ef)[[2]], -3)
  expect_true(all(g@bread[2, ] == 0))
  expect_true(is_root(ef, theta))
  not_converged <- Filter(
    function(cnd) inherits(cnd, "deli_solver_not_converged"),
    seen
  )
  expect_length(not_converged, 1L)
  reported <- conditionMessage(not_converged[[1]])
  expect_match(
    reported,
    "Moment condition 2 does not move when any parameter does",
    fixed = TRUE
  )
  expect_match(reported, "it sums to -3 rather than to zero", fixed = TRUE)
  # The flat branch rather than a reading of the contributions, which cancel
  # well enough here to say nothing.
  expect_no_match(reported, "do not cancel", fixed = TRUE)
})

# A flat equation is judged by the size of its summed score against the size of
# the one it had at the starting values, and a size is all that comparison reads.
# The known miss mode is a trade across the root: an equation the solver moved
# from one side of its root to the other, ending as far from it as it began,
# reads as unchanged and is not reported. The comparison is what keeps a median
# equation quiet where the caller started it at the sample median and the fit
# held it there, which is the documented way to use one, so this is the price of
# that and is recorded here rather than fixed.

test_that("a flat equation traded across its root is not judged", {
  # Twelve observations of a two-equation stack, the second of them flat: its
  # row of the bread is identically zero, so nothing about it moves when a
  # parameter does. Nine contributions of 1 against three of -1 leave it summing
  # to 6.
  ef <- rbind(rep(0.5, 12), c(rep(1, 9), rep(-1, 3)))
  bread <- rbind(c(1, 0.5), c(0, 0))

  # Started at -6, which is the same distance from the root on the other side.
  # The trade is invisible to the comparison and the point passes.
  expect_null(flat_equation(ef, bread, init_score = c(0, -6)))

  # A score that grew rather than changed sides is judged, whichever sign it
  # carries at either end.
  found <- flat_equation(ef, bread, init_score = c(0, -5))
  expect_equal(found$row, 2L)
  expect_equal(found$score, 6)
  expect_true(found$flat)
})

# ---- a runaway estimate as the scale it is measured against -------------------
#
# The one-sided reading measures the mean contribution of an equation that
# cannot cancel against the scale of the problem, and that scale takes in the
# magnitude of the estimates. The term is there for a stack that is nothing but
# such an equation, where the estimates are the only scale available, and it is
# self-defeating exactly when the estimates are the thing that went wrong: a
# solver that runs a parameter away supplies a scale large enough to excuse
# whatever it left behind.

test_that("a runaway estimate does not supply the scale it is judged against", {
  # inverse_logit is bounded by one, so no value of theta brings this equation
  # near zero: it sums to -160 wherever the solver stops. multiroot runs theta
  # out to 3.4e8 and reports success there. The equation is one-sided, so it is
  # measured against the scale of the problem, and the runaway estimate is
  # itself the largest thing in that scale: 160 / 40 / 3.4e8 reads as 1.2e-8,
  # under the ceiling.
  psi <- function(theta) matrix(inverse_logit(theta[1]) - rep(5, 40), nrow = 1)
  seen <- collect_warnings({
    m <- estimate(MEstimator(stacked_equations = psi, init = 0))
  })
  expect_length(seen, 1L)
  expect_length(
    Filter(function(cnd) inherits(cnd, "deli_solver_not_converged"), seen),
    1L
  )
  theta <- unname(coef(m))
  expect_gt(abs(theta), 1e6)
  # Nothing else reaches it. A saturated logistic has a bread of zeros, so no
  # Newton step can be taken, and the score the solver left is smaller than the
  # one it was handed, so the flat-bread reading owes it the same excuse it owes
  # an equation left where the caller put it.
  expect_true(all(m@bread == 0))
  expect_equal(sum(psi(theta)), -160)
  expect_lt(abs(sum(psi(theta))), abs(sum(psi(0))))
})

test_that("the same runaway equation is judged on the GMM path", {
  # The asymmetry is the evidence that the scale rather than the equation is
  # what goes wrong above: minimizing the same stack leaves theta at 12, which
  # is too small to excuse anything, and the moment is reported.
  psi <- function(theta) matrix(inverse_logit(theta[1]) - rep(5, 40), nrow = 1)
  seen <- collect_warnings(
    estimate(GMMEstimator(stacked_equations = psi, init = 0))
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  expect_match(
    conditionMessage(seen[[1]]),
    "does not vanish, leaving a summed moment of -160",
    fixed = TRUE
  )
})

test_that("a solved relation is still measured against the scale of the estimates", {
  # The guard the term is owed. A mean equation on data of magnitude 4e8 and a
  # contrast demanding the two parameters agree, handed to a solver whose
  # convergence test is relative to the estimates, which is what every parameter
  # tolerance on offer is written in terms of. It stops where the two agree to
  # ten significant digits, leaving the contrast summing to 8 rather than to
  # zero while every contribution in the stack is of order one.
  set.seed(9)
  n <- 200
  y <- 4e8 + stats::rnorm(n)
  psi <- function(theta) {
    rbind(
      y - theta[1],
      rep(theta[1] - theta[2], n)
    )
  }
  relative_tolerance <- function(stacked_equations, init) {
    root <- mean(y)
    c(root, root * (1 - 1e-10))
  }
  expect_no_warning({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(4e8, 4e8)),
      solver = relative_tolerance
    )
  })
  theta <- unname(coef(m))
  ef <- psi(theta)
  expect_true(is_root(ef, theta))
  # The contrast is one-sided and its score is well above the floor, so the
  # reading is made rather than excused.
  expect_equal(sum(ef[2, ] < 0), 0L)
  expect_equal(rowSums(ef)[[2]], 8, tolerance = 1e-6)
  expect_gt(abs(rowSums(ef)[[2]]), score_floor)
  # Measured against the contributions alone the contrast is out by a factor of
  # 144, so a scale carrying no reading of the parameters at all turns this
  # solved fit into a false alarm.
  expect_gt(
    abs(rowSums(ef)[[2]]) / n / (one_sided_ceiling * max(1, abs(ef))),
    100
  )
})

# ---- estimates that run away while the score does not ------------------------
#
# A flat estimating function drives the summed score below any floor while the
# parameters run away rather than settle, so every reading written in terms of
# the size of the residual is met at a point nobody would accept. Catching it
# needs a reading of the estimates themselves, and the guards below are what
# keeps such a reading from reporting a fit that legitimately traveled a long
# way to a root.

test_that("a separated logistic that runs away is not reported as solved", {
  # A covariate that separates the outcome exactly. The likelihood has no
  # maximum: the slope grows without bound and the score falls with it, so
  # multiroot's own atol test is met at a slope of 11007 with a score of 7.6e-06.
  set.seed(42)
  n <- 200
  x <- stats::rnorm(n)
  y <- as.numeric(x > 0)
  design <- cbind(1, x)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = y, model = "logistic")
  }
  seen <- collect_warnings({
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0)))
  })
  expect_length(seen, 1L)
  expect_length(
    Filter(function(cnd) inherits(cnd, "deli_solver_not_converged"), seen),
    1L
  )
  theta <- unname(coef(m))
  # Where a runaway stops is a floating-point accident, so the slope is bounded
  # rather than pinned. It reads 11007 on macOS and in both Linux containers,
  # which agree to six digits, and the bound sits an order of magnitude under
  # that while staying far above any slope this data could support.
  expect_gt(theta[[2]], 1e3)
  # No reading of the contributions sees it: the score is under the floor, so
  # they say nothing at all, and the solver reported success rather than a
  # convergence test that failed. What is left is the Newton step, and it is
  # bounded rather than pinned because the bread is near singular where the
  # parameters run away, so its inverse magnifies the last digits of the
  # estimates. The step reads 0.102 here and in both Linux containers and 0.091
  # on the R 4.4 Linux CI runner, and perturbing the estimates in their eighth
  # digit moves it over that same range. The bound is an order of magnitude under
  # the smallest of those readings, an order above the ceiling the step is judged
  # against, and six orders above the fit that travels a long way to a root
  # below, which is the contrast this reading exists to draw.
  ef <- psi(theta)
  expect_lt(max(abs(rowSums(ef))), score_floor)
  expect_true(is_root(ef, theta))
  expect_gt(
    relative_newton_step(unname(m@bread), rowSums(ef) / n, theta),
    0.01
  )
})

test_that("a logistic fit on overlapping data stays quiet", {
  set.seed(7)
  n <- 200
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(-0.3 + 1.1 * x))
  design <- cbind(1, x)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = y, model = "logistic")
  }
  expect_no_warning({
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0)))
  })
  expect_equal(
    unname(coef(m)),
    unname(stats::coef(stats::glm(y ~ x, family = stats::binomial()))),
    tolerance = 1e-6
  )
})

test_that("a fit that travels a long way to a root stays quiet", {
  # The estimate moves a million units from the starting value and lands on the
  # sample mean, so a bound on how far a fit may travel from where it started
  # would report this one. Where it stopped is what separates it from the
  # separated logistic above: the Newton step there is a tenth, and here it is
  # 2.4e-09 on macOS and 6.2e-12 in the Linux containers, so the bound below
  # clears the larger of them by better than two orders of magnitude.
  set.seed(12)
  n <- 200
  y <- 1e6 + stats::rnorm(n)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  expect_no_warning({
    m <- estimate(MEstimator(stacked_equations = psi, init = 0))
  })
  theta <- unname(coef(m))
  expect_equal(theta, mean(y), tolerance = 1e-6)
  expect_gt(abs(theta), 1e5)
  expect_lt(
    relative_newton_step(unname(m@bread), rowSums(psi(theta)) / n, theta),
    1e-6
  )
})

# ---- a design that does not identify the parameters --------------------------
#
# A design whose columns are linearly dependent leaves the estimating equations
# with no unique root: the returned point solves them, and so does every point
# along the direction the design cannot see. What a solver makes of such a system
# is not fixed. Faced with a singular Jacobian it may walk to a root and report
# nothing of its own, or stop where it started and report a convergence test that
# failed, and which of the two happens depends on the LAPACK the platform's R was
# built against. A report about the search is true of the search and beside the
# point about the problem on either branch, so what the tests below hold to is
# the wording rather than the branch: however the solve went, the fit names the
# design as rank deficient and the parameters as not identified. The tests that
# read anything more than the wording pin one branch deliberately, from starting
# values that already solve the equations or through a solver that reports a
# convergence test it met.

# A linear regression whose design repeats its second column, so the last two
# coefficients are identified only by their sum.
duplicated_column_psi <- function() {
  set.seed(3)
  n <- 200
  x1 <- stats::rnorm(n)
  y <- 1 + 2 * x1 + stats::rnorm(n)
  design <- cbind(1, x1, x1)
  function(theta) ee_regression(theta, X = design, y = y, model = "linear")
}

# The same data and shape with a third column of its own.
full_rank_psi <- function() {
  set.seed(3)
  n <- 200
  x1 <- stats::rnorm(n)
  y <- 1 + 2 * x1 + stats::rnorm(n)
  x2 <- stats::rnorm(n)
  design <- cbind(1, x1, x2)
  function(theta) ee_regression(theta, X = design, y = y, model = "linear")
}

# The reading itself, at the shapes a bread can arrive in. Each guard is there
# for a bread some fit in this suite produces, and a guard that stopped holding
# would turn one of those into a false alarm rather than into an error, so each
# is pinned on its own.

test_that("rank_deficient() reads only a square, finite bread", {
  expect_false(rank_deficient(rbind(c(1, 0), c(0, 1))))
  expect_true(rank_deficient(rbind(c(1, 2), c(2, 4))))
  # A subset fit hands over the block of the bread its solver worked on, and a
  # caller that passed the whole of it would be asking about a shape this
  # reading has no answer for.
  expect_false(rank_deficient(rbind(c(1, 2, 3), c(2, 4, 6))))
  expect_false(rank_deficient(matrix(numeric(0), nrow = 0, ncol = 0)))
  # qr() cannot factor a matrix carrying values that are not finite, and a
  # derivative that is not finite is a different failure with a different
  # remedy.
  expect_false(rank_deficient(rbind(c(NA_real_, 1), c(1, 1))))
  expect_false(rank_deficient(rbind(c(Inf, 1), c(1, 1))))
  # A bread of zeros is rank deficient in every direction there is.
  expect_true(rank_deficient(matrix(0, nrow = 2, ncol = 2)))
})

test_that("not_identified() leaves a bread that lost its rank to a flat row", {
  # The duplicated-column shape: dependent rows, none of them zero.
  expect_true(not_identified(rbind(c(1, 2), c(2, 4))))
  # The non-differentiable shape: one equation that does not move at all. The
  # remedy this reading names is not the remedy for that, so it says nothing.
  expect_false(not_identified(rbind(c(1, 1), c(0, 0))))
  expect_false(not_identified(matrix(0, nrow = 2, ncol = 2)))
  expect_false(not_identified(rbind(c(1, 0), c(0, 1))))
})

test_that("a rank-deficient design is reported as such", {
  psi <- duplicated_column_psi()
  # Nothing here reads the point the fit returned. Which point that is, and
  # whether the solver was satisfied with it, is the platform's business: from
  # these starting values the search either reaches a root of the singular
  # system or gives up at its Jacobian and hands the starting values back. One
  # warning comes out of the fit on either branch, and it is the one that names
  # the design.
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = c(0, 0, 0)))
  )
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  reported <- conditionMessage(seen[[1]])
  expect_match(reported, "rank deficient", fixed = TRUE)
  expect_match(reported, "not identified", fixed = TRUE)
})

test_that("the rank-deficient report does not read as a search that stopped short", {
  psi <- duplicated_column_psi()
  # Starting values that already solve the equations, so the solver accepts them
  # where it finds them and states nothing about a search on any platform. What
  # is left to word is the reading of the returned point, and a search that
  # stopped short is not what that reading found.
  init <- c(1.044995, 0, 1.935798)
  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = init))
  )
  expect_length(seen, 1L)
  expect_no_match(
    conditionMessage(seen[[1]]),
    "did not converge to a root of the estimating equations",
    fixed = TRUE
  )
})

test_that("a rank-deficient design is reported from starting values that solve it", {
  psi <- duplicated_column_psi()
  # The values the fit above returns, handed back to it as starting values. The
  # solver accepts them where it finds them and reports nothing at all, so
  # nothing about the search is left to word: the report has to come from the
  # reading of the returned point.
  init <- c(1.044995, 0, 1.935798)
  seen <- collect_warnings({
    m <- estimate(MEstimator(stacked_equations = psi, init = init))
  })
  expect_length(seen, 1L)
  reported <- conditionMessage(seen[[1]])
  expect_match(reported, "rank deficient", fixed = TRUE)
  expect_match(reported, "not identified", fixed = TRUE)
  expect_true(is_root(psi(unname(coef(m))), unname(coef(m))))
})

test_that("a rank-deficient design is reported whichever solver returned the point", {
  skip_if_not_installed("minpack.lm")
  psi <- duplicated_column_psi()
  # minpack.lm reports a convergence test that was met and says nothing more,
  # which is the same silence a custom solver returns. The design is no less
  # rank deficient for it.
  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = c(0, 0, 0)),
      solver = "lm"
    )
  })
  expect_length(seen, 1L)
  expect_match(conditionMessage(seen[[1]]), "not identified", fixed = TRUE)
  theta <- unname(coef(m))
  expect_equal(theta[[2]] + theta[[3]], 1.9358, tolerance = 1e-4)
})

test_that("a full-rank design of the same shape stays quiet", {
  psi <- full_rank_psi()
  expect_no_warning({
    m <- estimate(MEstimator(stacked_equations = psi, init = c(0, 0, 0)))
  })
  expect_equal(
    unname(coef(m)),
    c(1.04503043, 1.93571355, 0.00753299),
    tolerance = 1e-6
  )
})

# ---- conditions crossing the solver boundary ---------------------------------
#
# One fit evaluates its estimating function many times: at the starting values,
# inside the solver, at the point the solver returned, and once or twice per
# parameter while the bread is built. A warning the estimating function raises is
# the user's business whichever of those evaluations raised it, so no solver
# branch may swallow one. What is not the user's business is a solver's own
# account of itself, which `solve_equations()` reads and rewords, and which is
# therefore recognized by what it says and muffled where it is raised.

# A psi whose warning names the evaluation that raised it, so no two of its
# warnings share a message and the de-duplication scope collapses none of them.
# The count delivered is then the count raised. The counter lives in an
# environment so the test can read the total back after the fit.
new_evaluation_counter <- function() {
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  counter
}

counting_warn_psi <- function(y, counter) {
  function(theta) {
    counter$n <- counter$n + 1L
    cli::cli_warn("evaluation {counter$n}")
    matrix(y - theta[1], nrow = 1)
  }
}

# A mean whose estimating function is otherwise quiet, small enough that every
# solver reaches the solution in a handful of steps.
warning_psi_data <- function() {
  set.seed(1)
  stats::rnorm(40)
}

test_that("the default solver delivers every warning its psi raises", {
  counter <- new_evaluation_counter()
  psi <- counting_warn_psi(warning_psi_data(), counter)

  seen <- collect_warnings(m_estimate(psi, init = 0))

  # Four evaluations happen outside the solver: one at the starting values, one
  # at the point returned, and two for the bread of a one-parameter fit. More
  # than four in total is what makes the count evidence about the evaluations
  # the solver made itself.
  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("evaluation", seq_len(counter$n))
  )
})

test_that("the nleqslv solver delivers every warning its psi raises", {
  skip_if_not_installed("nleqslv")
  counter <- new_evaluation_counter()
  psi <- counting_warn_psi(warning_psi_data(), counter)

  seen <- collect_warnings(m_estimate(psi, init = 0, solver = "nleqslv"))

  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("evaluation", seq_len(counter$n))
  )
})

test_that("the lm solver delivers every warning its psi raises", {
  skip_if_not_installed("minpack.lm")
  counter <- new_evaluation_counter()
  psi <- counting_warn_psi(warning_psi_data(), counter)

  seen <- collect_warnings(m_estimate(psi, init = 0, solver = "lm"))

  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("evaluation", seq_len(counter$n))
  )
})

# multiroot reports its own convergence test failing as a warning reading
# "steady-state not reached", and reports a Jacobian it cannot factor as one
# reading "error during factorisation of matrix (dgefa)". Both are the solver
# talking about itself, both are read by `solve_equations()` and reworded, and
# neither may reach the user twice. Both are quoted from rootSolve, spelling
# included, so an en-US sweep that respelled "factorisation" here would be
# rewriting another package's message.
test_that("multiroot's own account of a diverging solve stays inside the solver", {
  psi <- beta_runaway_psi()

  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = c(0, 0, log(10))))
  )

  reported <- paste(
    vapply(seen, conditionMessage, character(1)),
    collapse = " "
  )
  expect_no_match(reported, "steady-state not reached", fixed = TRUE)
  expect_no_match(reported, "dgefa", fixed = TRUE)
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
})

test_that("multiroot's own account of an exhausted budget stays inside the solver", {
  psi <- ratio_psi(3)

  seen <- collect_warnings(
    estimate(
      MEstimator(stacked_equations = psi, init = c(1, 1, 1)),
      maxiter = 1
    )
  )

  reported <- paste(
    vapply(seen, conditionMessage, character(1)),
    collapse = " "
  )
  expect_no_match(reported, "steady-state not reached", fixed = TRUE)
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
})

# The de-duplication scope keys on the class vector and the message, so an
# estimating function that raises one wording at every evaluation reports once
# however many of those evaluations reach the caller. Delivering the solver's own
# evaluations therefore leaves the warn-once promise where it was.
test_that("one wording raised at every evaluation still reports once", {
  skip_if_not_installed("nleqslv")
  skip_if_not_installed("minpack.lm")
  y <- warning_psi_data()

  for (solver in c("rootSolve", "nleqslv", "lm")) {
    counter <- new_evaluation_counter()
    psi <- function(theta) {
      counter$n <- counter$n + 1L
      cli::cli_warn("the estimating equation is not differentiable")
      matrix(y - theta[1], nrow = 1)
    }

    seen <- collect_warnings(m_estimate(psi, init = 0, solver = solver))

    expect_gt(counter$n, 4L)
    expect_length(seen, 1L)
  }
})

# A solver's own report is recognized by what it says, because multiroot gives
# nothing else to recognize it by: its returned list carries no status flag. What
# says those words is not always the solver, though. An estimating function of
# the caller's own, or a solver of theirs reaching for rootSolve from inside one
# of its evaluations, reports in the same words from code the package's solve
# did not run. Muffling those costs the caller a report that is theirs and
# credits the solve with a failure it did not have, so what is muffled is decided
# by whose code was running as well as by what it said.

# The wording is numbered for the reason counting_warn_psi() numbers its own: no
# two of these warnings share a message, so the de-duplication scope collapses
# none of them and the count delivered is the count raised. Every one of them
# reads as the solver's own report, and every one of them is the caller's.
reporting_psi <- function(y, counter, report) {
  function(theta) {
    counter$n <- counter$n + 1L
    cli::cli_warn(paste(report, counter$n))
    matrix(y - theta[1], nrow = 1)
  }
}

test_that("a multiroot report the estimating function raised reaches the caller", {
  counter <- new_evaluation_counter()
  y <- warning_psi_data()
  psi <- reporting_psi(y, counter, "steady-state not reached")

  seen <- collect_warnings({
    m <- m_estimate(psi, init = 0)
  })

  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("steady-state not reached", seq_len(counter$n))
  )
  # And the solve itself is not reported as a failure on the strength of them.
  # It reached the mean, which is the whole of what it had to do.
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-8)
})

test_that("an inner solve inside an outer one still reads its own reports", {
  skip_if_not_installed("minpack.lm")
  set.seed(7)
  outer_y <- stats::rnorm(40, mean = 1)
  inner <- beta_runaway_psi()
  # The two-stage estimator the nested-solve section below describes: an inner
  # fit on the default solver, started from an estimating function an outer
  # solver that is not rootSolve is evaluating. The inner fit diverges, so it has
  # a report of its own to make, and it is the inner solve's own machinery that
  # has to make it.
  psi <- function(theta) {
    estimate(MEstimator(stacked_equations = inner, init = c(0, 0, log(10))))
    matrix(outer_y - theta[1], nrow = 1)
  }

  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = psi, init = 0),
      solver = "lm"
    )
  })

  # One warning, in deli's words, with multiroot's own account of itself left
  # inside the solve that produced it.
  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  reported <- conditionMessage(seen[[1]])
  expect_no_match(reported, "steady-state not reached", fixed = TRUE)
  expect_no_match(reported, "dgefa", fixed = TRUE)
  expect_equal(unname(coef(m)), mean(outer_y), tolerance = 1e-8)
})

test_that("an nls.lm report the estimating function raised reaches the caller", {
  skip_if_not_installed("minpack.lm")
  counter <- new_evaluation_counter()
  y <- warning_psi_data()
  psi <- reporting_psi(y, counter, "lmdif: info =")

  seen <- collect_warnings({
    m <- m_estimate(psi, init = 0, solver = "lm")
  })

  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("lmdif: info =", seq_len(counter$n))
  )
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-8)
})

# ---- nested solves -----------------------------------------------------------
#
# rootSolve::multiroot() cannot be called from inside itself. Its C code keeps
# the environment the estimating function is evaluated in in a single slot, so an
# inner call overwrites what the outer call is still using and the outer solve
# carries on over corrupted state. What that looks like depends on the shape of
# the two problems: the outer solve may fail out of the C code with a type error,
# or it may return quietly with the inner fit's estimates in place of its own.
# Neither may be allowed to happen, so a rootSolve solve started while another
# one is running is refused outright.

# Two independent samples with different means, so an outer fit that comes back
# with the inner fit's answer is visible in the estimate rather than only in the
# diagnostics.
nested_solve_data <- function() {
  set.seed(11)
  list(
    inner = stats::rnorm(40, mean = 2),
    outer = stats::rnorm(40, mean = 1),
    z = stats::rnorm(40)
  )
}

# An outer estimating function that fits an inner M-estimator before returning
# its own contributions. `inner_solver` names the solver the inner fit uses, and
# NULL leaves it on the default.
nested_mean_psi <- function(d, inner_solver = NULL) {
  inner <- function(theta) matrix(d$inner - theta[1], nrow = 1)
  function(theta) {
    estimate(
      MEstimator(stacked_equations = inner, init = 0),
      solver = inner_solver
    )
    matrix(d$outer - theta[1], nrow = 1)
  }
}

# The same nesting with a two-parameter outer problem, which corrupts the C
# state differently from the one-parameter one.
nested_regression_psi <- function(d) {
  inner <- function(theta) matrix(d$inner - theta[1], nrow = 1)
  x <- cbind(1, d$z)
  function(theta) {
    estimate(MEstimator(stacked_equations = inner, init = 0))
    t(x * as.vector(d$outer - x %*% theta))
  }
}

test_that("a rootSolve solve inside a rootSolve solve is refused", {
  d <- nested_solve_data()

  err <- expect_error(
    m_estimate(nested_mean_psi(d), init = 0),
    class = "deli_nested_solver_error"
  )

  reported <- gsub("\\s+", " ", conditionMessage(err))
  expect_match(reported, "rootSolve", fixed = TRUE)
  expect_match(reported, "nested", fixed = TRUE)
  # The remedy is a different solver for one of the two fits, and the message
  # has to name one that will work.
  expect_match(reported, "nleqslv", fixed = TRUE)
})

test_that("a nested rootSolve solve is refused whatever the outer problem", {
  d <- nested_solve_data()

  expect_error(
    m_estimate(nested_regression_psi(d), init = c(0, 0)),
    class = "deli_nested_solver_error"
  )
})

test_that("a plain fit still solves after a nested solve has been refused", {
  d <- nested_solve_data()
  expect_error(
    m_estimate(nested_mean_psi(d), init = 0),
    class = "deli_nested_solver_error"
  )

  m <- m_estimate(
    function(theta) matrix(d$outer - theta[1], nrow = 1),
    init = 0
  )

  expect_equal(unname(coef(m)), mean(d$outer), tolerance = 1e-8)
})

# Only a rootSolve solve started inside another one is refused. An inner
# rootSolve fit under an outer fit on any other solver is the pattern a
# two-stage estimator is written in, and it works.
test_that("a rootSolve solve inside a fit on another solver is allowed", {
  skip_if_not_installed("nleqslv")
  skip_if_not_installed("minpack.lm")
  d <- nested_solve_data()
  least_squares <- function(stacked_equations, init) {
    stats::optim(
      init,
      function(theta) sum(stacked_equations(theta)^2),
      method = "BFGS"
    )$par
  }

  for (solver in list(least_squares, "nleqslv", "lm")) {
    m <- estimate(
      MEstimator(stacked_equations = nested_mean_psi(d), init = 0),
      solver = solver
    )
    expect_equal(unname(coef(m)), mean(d$outer), tolerance = 1e-6)
  }
})

test_that("a fit on another solver inside a rootSolve solve is allowed", {
  skip_if_not_installed("nleqslv")
  skip_if_not_installed("minpack.lm")
  d <- nested_solve_data()

  for (solver in c("nleqslv", "lm")) {
    m <- m_estimate(nested_mean_psi(d, inner_solver = solver), init = 0)
    expect_equal(unname(coef(m)), mean(d$outer), tolerance = 1e-8)
  }
})

# ---- the lm solver's own warnings --------------------------------------------
#
# minpack.lm::nls.lm() reports an exhausted iteration budget twice: once as a
# bare warning out of lmdif carrying the MINPACK info code, and once in its
# return value, which is where `solve_equations()` reads it and words deli's own
# report. Only the second is the user's business, so the first is muffled by
# what it says, exactly as multiroot's own reports are.

test_that("a capped lm fit reports deli's warning alone", {
  skip_if_not_installed("minpack.lm")
  ref <- load_fixture("ee_solver_lm_logistic")
  psi <- function(theta) {
    ee_regression(theta, X = ref$X, y = ref$y, model = "logistic")
  }

  seen <- collect_warnings(
    estimate(
      MEstimator(stacked_equations = psi, init = ref$init),
      solver = "lm",
      maxiter = 1
    )
  )

  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  reported <- paste(
    vapply(seen, conditionMessage, character(1)),
    collapse = " "
  )
  expect_no_match(reported, "lmdif", fixed = TRUE)
  expect_no_match(reported, "info = -1", fixed = TRUE)
})

test_that("muffling lm's own warning leaves its psi's warnings alone", {
  skip_if_not_installed("minpack.lm")
  ref <- load_fixture("ee_solver_lm_logistic")
  counter <- new_evaluation_counter()
  psi <- function(theta) {
    counter$n <- counter$n + 1L
    cli::cli_warn("evaluation {counter$n}")
    ee_regression(theta, X = ref$X, y = ref$y, model = "logistic")
  }

  seen <- collect_warnings(
    estimate(
      MEstimator(stacked_equations = psi, init = ref$init),
      solver = "lm"
    )
  )

  expect_gt(counter$n, 4L)
  expect_equal(
    vapply(seen, conditionMessage, character(1)),
    paste("evaluation", seq_len(counter$n))
  )
})

# ---- an ill-conditioned nleqslv Jacobian -------------------------------------
#
# nleqslv reports termination code 5 where the Jacobian of the estimating
# equations is too ill-conditioned for it to take another step. That is not a
# report that the root was missed: a stack whose two blocks differ by orders of
# magnitude reaches code 5 at the exact solution, with the returned values equal
# to the ones the fit started from to the last digit. What an ill-conditioned
# Jacobian does put in doubt is the bread matrix built there, and so the
# variance, which is what the report has to say.

# Two independent means whose estimating functions differ by fifteen orders of
# magnitude. The scaling changes neither solution, and each block on its own is
# well behaved, so the only thing wrong at the solution is the conditioning.
ill_conditioned_means <- function() {
  set.seed(5)
  list(
    first = stats::rnorm(40, mean = 3),
    second = stats::rnorm(40, mean = 7)
  )
}

ill_conditioned_psi <- function(d) {
  function(theta) {
    rbind(1e9 * (d$first - theta[1]), 1e-6 * (d$second - theta[2]))
  }
}

test_that("an ill-conditioned nleqslv Jacobian is reported as a doubt about the variance", {
  skip_if_not_installed("nleqslv")
  d <- ill_conditioned_means()
  solution <- c(mean(d$first), mean(d$second))

  seen <- collect_warnings({
    m <- estimate(
      MEstimator(stacked_equations = ill_conditioned_psi(d), init = solution),
      solver = "nleqslv"
    )
  })

  expect_length(seen, 1L)
  expect_s3_class(seen[[1]], "deli_solver_not_converged")
  reported <- gsub("\\s+", " ", conditionMessage(seen[[1]]))
  expect_match(reported, "ill-conditioned", fixed = TRUE)
  expect_match(reported, "variance", fixed = TRUE)
  # The point being reported on is the solution the fit started from, unmoved,
  # so nothing in the report may read as a search that missed it.
  expect_equal(unname(m@theta), solution, tolerance = 1e-15)
  expect_no_match(reported, "did not converge", fixed = TRUE)
  expect_no_match(reported, "stopped without solving", fixed = TRUE)
  # rootSolve is the solver that returns a spurious root silently, so no
  # non-convergence message may send a user to it as the remedy.
  expect_no_match(reported, "rootSolve", fixed = TRUE)
})

test_that("a stalled nleqslv fit keeps its own wording", {
  skip_if_not_installed("nleqslv")
  # theta^2 + 1 has no real root, and nleqslv stalls rather than failing hard.
  psi <- function(theta) matrix(rep((theta[1]^2 + 1) / 10, 10), nrow = 1)

  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = 1), solver = "nleqslv")
  )

  expect_length(seen, 1L)
  reported <- gsub("\\s+", " ", conditionMessage(seen[[1]]))
  expect_match(
    reported,
    "stopped without solving the estimating equations (code 3)",
    fixed = TRUE
  )
  expect_no_match(reported, "ill-conditioned", fixed = TRUE)
})

test_that("a singular nleqslv Jacobian keeps the non-convergence wording", {
  skip_if_not_installed("nleqslv")
  y <- rep(5, 40)
  psi <- function(theta) matrix(inverse_logit(theta[1]) - y, nrow = 1)

  seen <- collect_warnings(
    estimate(MEstimator(stacked_equations = psi, init = 0), solver = "nleqslv")
  )

  expect_length(seen, 1L)
  reported <- gsub("\\s+", " ", conditionMessage(seen[[1]]))
  expect_match(reported, "nleqslv did not converge (code 6)", fixed = TRUE)
  expect_no_match(reported, "ill-conditioned", fixed = TRUE)
})

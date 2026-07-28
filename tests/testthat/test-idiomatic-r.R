# Tests for idiomatic R improvements: parameter naming, S3 generics,
# m_estimate/gmm_estimate wrappers, and broom tidiers

# ---- Helpers -----------------------------------------------------------------

make_fitted_mean <- function(named = FALSE) {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  init <- if (named) c(mean = 0) else c(0)
  m <- MEstimator(stacked_equations = psi, init = init)
  estimate(m)
}

make_fitted_regression <- function() {
  y <- mtcars$mpg
  X <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m <- MEstimator(
    stacked_equations = psi,
    init = c(intercept = 0, wt = 0, hp = 0)
  )
  estimate(m)
}

# ---- Phase 1: Parameter naming -----------------------------------------------

test_that("named init propagates to theta", {
  m <- make_fitted_mean(named = TRUE)
  expect_equal(names(m@theta), "mean")
})

test_that("unnamed init gets auto-generated names", {
  m <- make_fitted_mean(named = FALSE)
  expect_equal(names(m@theta), "theta_1")
})

test_that("named init propagates to variance dimnames", {
  m <- make_fitted_mean(named = TRUE)
  expect_equal(dimnames(m@variance), list("mean", "mean"))
})

test_that("multi-parameter named init propagates correctly", {
  m <- make_fitted_regression()
  expect_equal(names(m@theta), c("intercept", "wt", "hp"))
  expect_equal(rownames(m@variance), c("intercept", "wt", "hp"))
  expect_equal(colnames(m@variance), c("intercept", "wt", "hp"))
})

test_that("names propagate to confidence_intervals", {
  m <- make_fitted_regression()
  ci <- confidence_intervals(m)
  expect_equal(rownames(ci), c("intercept", "wt", "hp"))
})

test_that("names propagate to z_scores and p_values", {
  m <- make_fitted_regression()
  z <- z_scores(m)
  p <- p_values(m)
  expect_equal(names(z), c("intercept", "wt", "hp"))
  expect_equal(names(p), c("intercept", "wt", "hp"))
})

test_that("constructor preserves names on init", {
  m <- MEstimator(
    stacked_equations = function(theta) matrix(1 - theta[1], nrow = 1),
    init = c(mu = 0)
  )
  expect_equal(names(m@init), "mu")
})

# ---- Phase 2: S3 generics ---------------------------------------------------

test_that("coef() returns named theta", {
  m <- make_fitted_regression()
  expect_equal(coef(m), m@theta)
  expect_equal(names(coef(m)), c("intercept", "wt", "hp"))
})

test_that("vcov() returns named variance matrix", {
  m <- make_fitted_regression()
  expect_equal(vcov(m), m@variance)
  expect_equal(rownames(vcov(m)), c("intercept", "wt", "hp"))
})

test_that("confint() works with default level", {
  m <- make_fitted_regression()
  ci <- confint(m)
  expect_equal(ncol(ci), 2)
  expect_equal(nrow(ci), 3)
  expect_equal(colnames(ci), c("lower", "upper"))
})

test_that("confint() works with custom level", {
  m <- make_fitted_regression()
  ci_99 <- confint(m, level = 0.99)
  ci_95 <- confint(m, level = 0.95)
  # 99% CI should be wider
  expect_true(
    ci_99[1, "upper"] - ci_99[1, "lower"] >
      ci_95[1, "upper"] - ci_95[1, "lower"]
  )
})

test_that("confint() subsets by parm name", {
  m <- make_fitted_regression()
  ci <- confint(m, parm = "wt")
  expect_equal(nrow(ci), 1)
  expect_equal(rownames(ci), "wt")
})

test_that("confint() subsets by parm index", {
  m <- make_fitted_regression()
  ci <- confint(m, parm = 1:2)
  expect_equal(nrow(ci), 2)
})

test_that("nobs() returns observation count", {
  m <- make_fitted_regression()
  expect_equal(nobs(m), 32L)
})

test_that("S3 generics error before estimation", {
  psi <- function(theta) matrix(1 - theta[1], nrow = 1)
  m <- MEstimator(stacked_equations = psi, init = c(0))
  expect_error(coef(m))
  expect_error(vcov(m))
  expect_error(confint(m))
  expect_error(nobs(m))
})

test_that("S3 generics work for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  g <- GMMEstimator(stacked_equations = psi, init = c(mean = 0))
  g <- estimate(g)
  expect_equal(names(coef(g)), "mean")
  expect_true(!is.null(vcov(g)))
  expect_equal(nobs(g), 5L)
  expect_equal(nrow(confint(g)), 1)
})

# ---- Phase 4: m_estimate() and gmm_estimate() --------------------------------

test_that("m_estimate() formula interface works", {
  m <- m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  expect_s3_class(m, "deli::MEstimator")
  expect_equal(names(coef(m)), c("(Intercept)", "wt", "hp"))
  expect_equal(nobs(m), 32L)
})

test_that("m_estimate() formula auto-generates named init", {
  m <- m_estimate(
    mpg ~ wt,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  expect_equal(names(coef(m)), c("(Intercept)", "wt"))
})

test_that("m_estimate() default interface works", {
  y <- c(1, 2, 3, 4, 5)
  m <- m_estimate(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(mean = 0)
  )
  expect_s3_class(m, "deli::MEstimator")
  expect_equal(names(coef(m)), "mean")
  expect_equal(unname(coef(m)), mean(y), tolerance = 1e-4)
})

test_that("m_estimate() formula matches manual MEstimator workflow", {
  # Formula interface
  m1 <- m_estimate(
    mpg ~ wt,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )

  # Manual workflow
  X <- model.matrix(mpg ~ wt, data = mtcars)
  y <- mtcars$mpg
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")
  m2 <- MEstimator(stacked_equations = psi, init = c(`(Intercept)` = 0, wt = 0))
  m2 <- estimate(m2)

  expect_equal(coef(m1), coef(m2), tolerance = 1e-6)
  expect_equal(vcov(m1), vcov(m2), tolerance = 1e-6)
})

test_that("m_estimate() NSE evaluates data columns", {
  skip_if_not_installed("MASS")
  # Create survival-like data for ee_aft
  set.seed(42)
  n <- 50
  d <- data.frame(
    time = rexp(n, rate = 0.5),
    x1 = rnorm(n),
    event = rbinom(n, 1, 0.7)
  )

  # NSE: event = event is evaluated in data context; the formula response
  # (time) is passed positionally and reaches ee_aft's `time` argument.
  # ee_aft with weibull needs ncol(X) + 1 params (extra log-scale)
  m <- m_estimate(
    time ~ x1,
    data = d,
    .ee = ee_aft,
    event = event,
    distribution = "weibull",
    init = c(`(Intercept)` = 0, x1 = 0, log_scale = 0)
  )
  expect_s3_class(m, "deli::MEstimator")
  expect_true(length(coef(m)) >= 2)
})

test_that("gmm_estimate() default interface works", {
  y <- c(1, 2, 3, 4, 5)
  g <- gmm_estimate(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(mean = 0)
  )
  expect_s3_class(g, "deli::GMMEstimator")
  expect_equal(names(coef(g)), "mean")
  expect_equal(unname(coef(g)), mean(y), tolerance = 1e-4)
})

# ---- subset and finite_correction forwarded through the formula interface ----
#
# The one-step wrappers create an MEstimator or GMMEstimator and hand the
# estimator constructor's subset and finite_correction options straight through.
# subset restricts root-finding to the named parameter indices, freezing the
# others at their init values (the variance estimator ignores subset).
# finite_correction rescales the meat matrix and switches confidence intervals
# and P-values to the t-distribution with df = n_obs - n_params.
#
# These tests pin that forwarding: each wrapper fit equals the equivalent
# manual estimate() fit carrying the same option, and each option changes the
# result relative to the fit without it. Both are formals of all four wrapper
# methods, so neither name reaches the estimating equation through `...`.

test_that("m_estimate() forwards finite_correction to t-distribution inference", {
  X <- model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")
  init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

  m_wrap <- m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear",
    finite_correction = "HC1"
  )
  m_manual <- estimate(MEstimator(psi, init = init, finite_correction = "HC1"))

  # Wrapper fit equals the manual estimator carrying the same correction.
  expect_equal(coef(m_wrap), coef(m_manual), tolerance = 1e-6)
  expect_equal(vcov(m_wrap), vcov(m_manual), tolerance = 1e-6)
  expect_equal(
    confidence_intervals(m_wrap),
    confidence_intervals(m_manual),
    tolerance = 1e-6
  )
  expect_equal(p_values(m_wrap), p_values(m_manual), tolerance = 1e-6)

  # The correction discriminates against the uncorrected (normal) fit: the
  # HC1 rescaling and heavier t-tails both widen every interval.
  m_unc <- m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  ci_corr <- confidence_intervals(m_wrap)
  ci_unc <- confidence_intervals(m_unc)
  expect_true(all(ci_corr[, "lower"] < ci_unc[, "lower"]))
  expect_true(all(ci_corr[, "upper"] > ci_unc[, "upper"]))
  expect_true(all(p_values(m_wrap) > p_values(m_unc)))
})

test_that("m_estimate() forwards subset to the root-finder", {
  X <- model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")
  init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

  m_wrap <- m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear",
    subset = 1L
  )
  m_manual <- estimate(MEstimator(psi, init = init, subset = 1L))

  # Wrapper fit equals the manual estimator solving the same parameter subset.
  expect_equal(coef(m_wrap), coef(m_manual), tolerance = 1e-6)
  expect_equal(vcov(m_wrap), vcov(m_manual), tolerance = 1e-6)

  # subset freezes wt and hp at their init values and solves only the
  # intercept, which becomes the marginal mean of mpg.
  expect_equal(unname(coef(m_wrap)[c("wt", "hp")]), c(0, 0))
  expect_equal(
    unname(coef(m_wrap)[["(Intercept)"]]),
    mean(mtcars$mpg),
    tolerance = 1e-4
  )

  # The subset fit differs materially from the full fit.
  m_full <- m_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  expect_true(
    abs(
      coef(m_wrap)[["(Intercept)"]] -
        coef(m_full)[["(Intercept)"]]
    ) >
      10
  )
})

test_that("gmm_estimate() forwards finite_correction to t-distribution inference", {
  X <- model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")
  init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

  g_wrap <- gmm_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear",
    finite_correction = "HC1"
  )
  g_manual <- estimate(GMMEstimator(
    psi,
    init = init,
    finite_correction = "HC1"
  ))

  expect_equal(coef(g_wrap), coef(g_manual), tolerance = 1e-5)
  expect_equal(vcov(g_wrap), vcov(g_manual), tolerance = 1e-5)
  expect_equal(
    confidence_intervals(g_wrap),
    confidence_intervals(g_manual),
    tolerance = 1e-5
  )
  expect_equal(p_values(g_wrap), p_values(g_manual), tolerance = 1e-5)

  g_unc <- gmm_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  ci_corr <- confidence_intervals(g_wrap)
  ci_unc <- confidence_intervals(g_unc)
  expect_true(all(ci_corr[, "lower"] < ci_unc[, "lower"]))
  expect_true(all(ci_corr[, "upper"] > ci_unc[, "upper"]))
  expect_true(all(p_values(g_wrap) > p_values(g_unc)))
})

test_that("gmm_estimate() forwards subset to the minimizer", {
  X <- model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  psi <- function(theta) ee_regression(theta, X = X, y = y, model = "linear")
  init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

  g_wrap <- gmm_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear",
    subset = 1L
  )
  g_manual <- estimate(GMMEstimator(psi, init = init, subset = 1L))

  expect_equal(coef(g_wrap), coef(g_manual), tolerance = 1e-5)
  expect_equal(vcov(g_wrap), vcov(g_manual), tolerance = 1e-5)

  # subset freezes wt and hp at their init values.
  expect_equal(unname(coef(g_wrap)[c("wt", "hp")]), c(0, 0))

  # The subset fit differs materially from the full fit.
  g_full <- gmm_estimate(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  expect_true(
    abs(
      coef(g_wrap)[["(Intercept)"]] -
        coef(g_full)[["(Intercept)"]]
    ) >
      10
  )
})

# ---- One-step wrappers reproduce the two-step form ---------------------------
#
# m_estimate() and gmm_estimate() build an estimator and call estimate() in a
# single call, handing every constructor argument and every estimate() argument
# straight through. The documentation leans on that equivalence: examples and
# articles use the one-step form everywhere the two-step form is not itself the
# subject, so an argument that stopped reaching the estimator would change
# published numbers without changing any visible code.
#
# Each grid entry carries two claims. The first is the equivalence itself,
# checked with expect_identical() rather than a tolerance so that a refactor
# quietly rerouting an argument is caught instead of absorbed by slack. The
# second is that the argument set moves the fit away from the same fixture's
# default. The second claim is what keeps the first from being vacuous: an
# argument the fixture cannot feel reproduces the default fit wherever the
# wrapper sends it, so the two forms agree whether or not the pass-through
# works. Arguments the linear mtcars fit cannot feel therefore get fixtures of
# their own further down: a nonlinear fit for the M estimator's maxiter and
# tolerance, a rank-deficient design for allow_pinv, and an over-identified
# system for the two over-identification controls.

# The arguments the estimator constructors take. Everything else in an argument
# set belongs to estimate().
constructor_arg_names <- c(
  "subset",
  "finite_correction",
  "overid_maxiter",
  "overid_tolerance"
)

fit_two_step <- function(constructor, psi, init, args) {
  constructor_args <- args[names(args) %in% constructor_arg_names]
  estimate_args <- args[!names(args) %in% constructor_arg_names]
  obj <- do.call(
    constructor,
    c(list(stacked_equations = psi, init = init), constructor_args)
  )
  do.call(estimate, c(list(obj), estimate_args))
}

same_fit <- function(a, b) {
  identical(coef(a), coef(b)) && identical(vcov(a), vcov(b))
}

expect_same_fit <- function(one_step, two_step, label) {
  expect_identical(coef(one_step), coef(two_step), info = label)
  expect_identical(vcov(one_step), vcov(two_step), info = label)
}

# One row of an argument grid: the arguments to pass, whether the fit they
# produce must differ from the fixture's default fit, and the warning they are
# expected to raise. An exhausted iteration budget is a non-convergence warning
# by construction, so the only way to observe maxiter is to catch that warning
# rather than compare bare fits.
arg_case <- function(args, differs = TRUE, warning = NULL) {
  list(args = args, differs = differs, warning = warning)
}

# Forces a fit, requiring the case's declared warning when it declares one.
# Anything undeclared propagates and fails the run.
fit_case <- function(fit, warning, label) {
  if (is.null(warning)) {
    return(fit)
  }
  # expect_warning() hands back the condition rather than the value, so the
  # fitted object is caught on the side.
  caught <- NULL
  expect_warning(caught <- fit, warning, info = label)
  caught
}

# Runs every case in a grid through both forms. `one_step` takes an argument
# list and returns the wrapper fit; `two_step` returns the build-then-estimate
# fit for the same arguments.
expect_grid_matches <- function(cases, one_step, two_step, default_fit) {
  for (label in names(cases)) {
    case <- cases[[label]]
    wrapped <- fit_case(one_step(case$args), case$warning, label)
    manual <- fit_case(two_step(case$args), case$warning, label)
    expect_same_fit(wrapped, manual, label)
    if (case$differs) {
      expect_false(same_fit(wrapped, default_fit), info = label)
    } else {
      expect_same_fit(wrapped, default_fit, label)
    }
  }
}

# Leave-one-out over a combined argument set. Dropping any member must change
# the fit; a member that survives its own removal is one the combined set never
# exercised, whatever the grid claims to cover.
expect_every_argument_felt <- function(fit_with, args) {
  full <- fit_with(args)
  for (name in names(args)) {
    expect_false(
      same_fit(fit_with(args[setdiff(names(args), name)]), full),
      info = name
    )
  }
}

# Linear regression on mtcars: three parameters, solved reliably from zero
# starting values by every solver under test. Being linear, it reaches its root
# in a couple of Newton steps far inside any iteration budget and below any
# tolerance, so maxiter and tolerance move to the logistic fixture below.
mtcars_regression_psi <- function() {
  X <- stats::model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  function(theta) ee_regression(theta, X = X, y = y, model = "linear")
}

mtcars_regression_init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

# Both combined sets below hold deriv_method at a finite-difference method,
# because "exact" ignores dx and would leave the dx member inert, and neither
# carries maxiter, because an exhausted budget stops the solver before
# tolerance can matter.
m_argument_sets <- list(
  "defaults" = arg_case(list(), differs = FALSE),
  "solver = NULL" = arg_case(list(solver = NULL), differs = FALSE),
  "solver = rootSolve" = arg_case(list(solver = "rootSolve"), differs = FALSE),
  "solver = lm" = arg_case(list(solver = "lm")),
  "deriv_method = capprox" = arg_case(
    list(deriv_method = "capprox"),
    differs = FALSE
  ),
  "deriv_method = exact" = arg_case(list(deriv_method = "exact")),
  "dx" = arg_case(list(dx = 1e-6)),
  "subset" = arg_case(list(subset = 1L)),
  "finite_correction = HC1" = arg_case(list(finite_correction = "HC1")),
  "every argument at once" = arg_case(list(
    solver = "lm",
    deriv_method = "fapprox",
    dx = 1e-7,
    subset = c(1L, 2L),
    finite_correction = "HC1"
  ))
)

# GMM minimizes rather than root-finds, so the solver names are stats::optim
# methods. "Nelder-Mead" is a poor minimizer for this badly scaled objective
# and wanders far from the least-squares answer, which is the point: it is a
# deterministic fit that the default "BFGS" cannot reproduce, so a solver that
# stopped reaching estimate() would show up immediately. Wandering that far is
# also a failed solve, so the case declares the moment-quality warning; the
# alternative would be picking a solver that happens not to trip it, which would
# leave the grid testing less than it does now.
# Minimizing is iterative even for a linear system, so unlike the M grid the
# GMM grid can hold maxiter and tolerance here.
gmm_argument_sets <- list(
  "defaults" = arg_case(list(), differs = FALSE),
  "solver = NULL" = arg_case(list(solver = NULL), differs = FALSE),
  "solver = BFGS" = arg_case(list(solver = "BFGS"), differs = FALSE),
  "solver = Nelder-Mead" = arg_case(
    list(solver = "Nelder-Mead"),
    warning = "not solved at the estimated values"
  ),
  "deriv_method = capprox" = arg_case(
    list(deriv_method = "capprox"),
    differs = FALSE
  ),
  "deriv_method = exact" = arg_case(list(deriv_method = "exact")),
  "maxiter" = arg_case(list(maxiter = 2L), warning = "did not converge"),
  "tolerance" = arg_case(list(tolerance = 1e-6)),
  "dx" = arg_case(list(dx = 1e-6)),
  "subset" = arg_case(list(subset = 1L)),
  "finite_correction = HC1" = arg_case(list(finite_correction = "HC1")),
  "every argument at once" = arg_case(list(
    solver = "Nelder-Mead",
    tolerance = 1e-6,
    deriv_method = "fapprox",
    dx = 1e-7,
    subset = c(1L, 2L),
    finite_correction = "HC1"
  ))
)

m_one_step_function <- function(psi, init) {
  function(args) {
    do.call(m_estimate, c(list(stacked_equations = psi, init = init), args))
  }
}

m_one_step_formula <- function(formula, data, ...) {
  fixed <- c(list(formula, data = data), list(...))
  function(args) do.call(m_estimate, c(fixed, args))
}

test_that("m_estimate() function interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()
  one_step <- m_one_step_function(psi, mtcars_regression_init)
  two_step <- function(args) {
    fit_two_step(MEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(m_argument_sets, one_step, two_step, one_step(list()))
})

test_that("m_estimate() formula interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()
  one_step <- m_one_step_formula(
    mpg ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "linear"
  )
  two_step <- function(args) {
    fit_two_step(MEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(m_argument_sets, one_step, two_step, one_step(list()))
})

test_that("gmm_estimate() function interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()
  one_step <- function(args) {
    do.call(
      gmm_estimate,
      c(list(stacked_equations = psi, init = mtcars_regression_init), args)
    )
  }
  two_step <- function(args) {
    fit_two_step(GMMEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(gmm_argument_sets, one_step, two_step, one_step(list()))
})

test_that("gmm_estimate() formula interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()
  one_step <- function(args) {
    do.call(
      gmm_estimate,
      c(
        list(
          mpg ~ wt + hp,
          data = mtcars,
          .ee = ee_regression,
          model = "linear"
        ),
        args
      )
    )
  }
  two_step <- function(args) {
    fit_two_step(GMMEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(gmm_argument_sets, one_step, two_step, one_step(list()))
})

test_that("the combined M argument set exercises every argument it carries", {
  psi <- mtcars_regression_psi()
  expect_every_argument_felt(
    m_one_step_function(psi, mtcars_regression_init),
    m_argument_sets[["every argument at once"]]$args
  )
})

test_that("the combined GMM argument set exercises every argument it carries", {
  psi <- mtcars_regression_psi()
  # The combined set holds "Nelder-Mead", which wanders far enough from the
  # least-squares answer to be a failed solve, and dropping "subset" puts that
  # fit back in front of the moment-quality check. This test compares fits, not
  # conditions; the warning itself is asserted in test-estimate-solvers.R.
  fit_with <- function(args) {
    suppressWarnings(
      do.call(
        gmm_estimate,
        c(list(stacked_equations = psi, init = mtcars_regression_init), args)
      ),
      classes = "deli_solver_not_converged"
    )
  }
  expect_every_argument_felt(
    fit_with,
    gmm_argument_sets[["every argument at once"]]$args
  )
})

# ---- maxiter and tolerance need a nonlinear fit ------------------------------
#
# Logistic regression on the same three-parameter shape gives rootSolve a
# genuine Newton sequence to walk, which is the only setting where an M
# estimator can feel its iteration budget or its convergence tolerance. A
# capped budget stops the solver short of the root and warns; a loose tolerance
# accepts an earlier iterate and stays silent, because solve_equations() warns
# on the budget rather than on the tolerance. The design carries the same wt
# and hp columns as the linear fixture, so it takes the same init vector.
mtcars_logistic_psi <- function() {
  X <- stats::model.matrix(am ~ wt + hp, data = mtcars)
  y <- mtcars$am
  function(theta) ee_regression(theta, X = X, y = y, model = "logistic")
}

m_iteration_argument_sets <- list(
  "defaults" = arg_case(list(), differs = FALSE),
  "maxiter" = arg_case(list(maxiter = 2L), warning = "did not converge"),
  "tolerance" = arg_case(list(tolerance = 1e-1))
)

test_that("m_estimate() function interface forwards maxiter and tolerance", {
  psi <- mtcars_logistic_psi()
  one_step <- m_one_step_function(psi, mtcars_regression_init)
  two_step <- function(args) {
    fit_two_step(MEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(
    m_iteration_argument_sets,
    one_step,
    two_step,
    one_step(list())
  )
})

test_that("m_estimate() formula interface forwards maxiter and tolerance", {
  psi <- mtcars_logistic_psi()
  one_step <- m_one_step_formula(
    am ~ wt + hp,
    data = mtcars,
    .ee = ee_regression,
    model = "logistic"
  )
  two_step <- function(args) {
    fit_two_step(MEstimator, psi, mtcars_regression_init, args)
  }

  expect_grid_matches(
    m_iteration_argument_sets,
    one_step,
    two_step,
    one_step(list())
  )
})

# ---- allow_pinv needs a singular bread ---------------------------------------
#
# allow_pinv chooses between solve() and MASS::ginv() in build_sandwich(), so it
# only does anything when the bread is singular. Duplicating a predictor makes
# the three-parameter bread rank two. Root-finding stays well posed because
# subset solves only the intercept and the wt slope, leaving the duplicate
# frozen at its zero start, which is the ordinary least-squares fit of mpg on
# wt; the bread is still assembled over all three parameters and so is
# singular. That the allow_pinv = FALSE fits error is itself the proof the
# fixture is doing its job.
collinear_regression_data <- function() {
  d <- mtcars[c("mpg", "wt")]
  d$wt_copy <- d$wt
  d
}

collinear_regression_psi <- function(d) {
  X <- stats::model.matrix(mpg ~ wt + wt_copy, data = d)
  y <- d$mpg
  function(theta) ee_regression(theta, X = X, y = y, model = "linear")
}

collinear_regression_init <- c(`(Intercept)` = 0, wt = 0, wt_copy = 0)

test_that("every wrapper interface forwards allow_pinv to the bread", {
  skip_if_not_installed("MASS")
  d <- collinear_regression_data()
  psi <- collinear_regression_psi(d)
  init <- collinear_regression_init
  # Solving only the identified parameters keeps the root-finder quiet while
  # the bread stays rank deficient.
  fixed <- list(subset = c(1L, 2L))
  formula_args <- list(
    mpg ~ wt + wt_copy,
    data = d,
    .ee = ee_regression,
    model = "linear"
  )

  wrappers <- list(
    "m_estimate() function" = list(
      one_step = function(args) {
        do.call(
          m_estimate,
          c(list(stacked_equations = psi, init = init), fixed, args)
        )
      },
      two_step = function(args) {
        fit_two_step(MEstimator, psi, init, c(fixed, args))
      }
    ),
    "m_estimate() formula" = list(
      one_step = function(args) {
        do.call(m_estimate, c(formula_args, fixed, args))
      },
      two_step = function(args) {
        fit_two_step(MEstimator, psi, init, c(fixed, args))
      }
    ),
    "gmm_estimate() function" = list(
      one_step = function(args) {
        do.call(
          gmm_estimate,
          c(list(stacked_equations = psi, init = init), fixed, args)
        )
      },
      two_step = function(args) {
        fit_two_step(GMMEstimator, psi, init, c(fixed, args))
      }
    ),
    "gmm_estimate() formula" = list(
      one_step = function(args) {
        do.call(gmm_estimate, c(formula_args, fixed, args))
      },
      two_step = function(args) {
        fit_two_step(GMMEstimator, psi, init, c(fixed, args))
      }
    )
  )

  for (label in names(wrappers)) {
    wrapper <- wrappers[[label]]
    # TRUE falls through to the pseudo-inverse and both forms agree.
    expect_same_fit(
      wrapper$one_step(list(allow_pinv = TRUE)),
      wrapper$two_step(list(allow_pinv = TRUE)),
      label
    )
    # FALSE leaves solve() to fail, in both forms alike.
    expect_error(
      wrapper$one_step(list(allow_pinv = FALSE)),
      "singular",
      info = label
    )
    expect_error(
      wrapper$two_step(list(allow_pinv = FALSE)),
      "singular",
      info = label
    )
  }
})

# ---- Over-identified GMM -----------------------------------------------------
#
# GMM earns its keep when there are more moment conditions than parameters. Two
# instruments for a single treatment effect give two estimating equations and
# one parameter, which MEstimator() rejects outright. The two-step weight matrix
# needs more than the default ten updating iterations to settle here, which also
# makes this the one system where overid_maxiter and overid_tolerance are
# observable in the fit.
#
# The moment conditions are written against a design matrix and a response so
# that the same function serves as the .ee of a formula fit.
# prepare_formula_psi() passes allow_over_identification = TRUE on the GMM path,
# so a .ee returning more rows than the design has columns reaches the updating
# loop from a formula as readily as from a closure.

make_iv_data <- function() {
  set.seed(42)
  n <- 200
  z1 <- stats::rbinom(n, 1, 0.5)
  z2 <- stats::rnorm(n)
  u <- stats::rnorm(n)
  a <- 0.5 * z1 + 0.3 * z2 + u + stats::rnorm(n)
  y <- 2 * a - u + stats::rnorm(n)
  data.frame(z1 = z1, z2 = z2, u = u, a = a, y = y)
}

iv_moment_conditions <- function(theta, X, y, z1, z2) {
  residual <- as.numeric(y - X %*% theta)
  rbind(z1 * residual, z2 * residual)
}

make_iv_psi <- function(d) {
  X <- stats::model.matrix(y ~ a - 1, data = d)
  function(theta) {
    iv_moment_conditions(theta, X = X, y = d$y, z1 = d$z1, z2 = d$z2)
  }
}

test_that("the IV moment conditions are over-identified", {
  d <- make_iv_data()
  psi <- make_iv_psi(d)

  # Two equations, one parameter. M-estimation has no root to find.
  expect_equal(dim(psi(c(effect = 0))), c(2L, 200L))
  expect_error(
    estimate(MEstimator(stacked_equations = psi, init = c(effect = 0))),
    regexp = "one estimating equation per parameter"
  )
})

test_that("gmm_estimate() solves the over-identified IV system", {
  d <- make_iv_data()
  psi <- make_iv_psi(d)

  g <- gmm_estimate(
    stacked_equations = psi,
    init = c(effect = 0),
    overid_maxiter = 200L
  )

  expect_s3_class(g, "deli::GMMEstimator")
  expect_equal(nobs(g), 200L)
  expect_equal(names(coef(g)), "effect")

  # The updating loop converged, so the weight matrix moved off the identity it
  # starts from.
  expect_false(isTRUE(all.equal(g@weight_matrix, diag(2))))

  # The instruments recover the treatment effect of 2 used to generate the
  # data, and the interval covers it.
  expect_equal(unname(coef(g)), 2, tolerance = 0.15)
  ci <- confint(g)
  expect_lt(ci[1, "lower"], 2)
  expect_gt(ci[1, "upper"], 2)

  # The confounded least-squares estimate is further from the truth, so the
  # instruments are doing the work rather than the data being easy.
  ols <- unname(coef(stats::lm(y ~ a - 1, data = d)))
  expect_lt(abs(unname(coef(g)) - 2), abs(ols - 2))
})

# The fully converged fit is the baseline the other two are measured against:
# zero updating iterations leaves the identity weight matrix in place, and a
# tolerance of 1 stops the updating after a single step. The package default of
# ten iterations does not settle this system and warns, so it cannot serve as
# the baseline.
overid_argument_sets <- list(
  "converged updating" = arg_case(
    list(overid_maxiter = 200L),
    differs = FALSE
  ),
  "no updating step" = arg_case(list(overid_maxiter = 0L)),
  "loose overid_tolerance" = arg_case(list(
    overid_maxiter = 200L,
    overid_tolerance = 1
  ))
)

test_that("gmm_estimate() over-identified fit reproduces the two-step fit", {
  d <- make_iv_data()
  psi <- make_iv_psi(d)
  init <- c(effect = 0)

  one_step <- function(args) {
    do.call(gmm_estimate, c(list(stacked_equations = psi, init = init), args))
  }
  two_step <- function(args) fit_two_step(GMMEstimator, psi, init, args)

  expect_grid_matches(
    overid_argument_sets,
    one_step,
    two_step,
    one_step(list(overid_maxiter = 200L))
  )
  expect_identical(
    one_step(list(overid_maxiter = 0L))@weight_matrix,
    diag(2)
  )
})

test_that("gmm_estimate() formula interface reaches the updating loop", {
  d <- make_iv_data()
  psi <- make_iv_psi(d)
  init <- c(a = 0)

  one_step <- function(args) {
    do.call(
      gmm_estimate,
      c(
        list(
          y ~ a - 1,
          data = d,
          .ee = iv_moment_conditions,
          z1 = d$z1,
          z2 = d$z2
        ),
        args
      )
    )
  }
  two_step <- function(args) fit_two_step(GMMEstimator, psi, init, args)

  # Two moment conditions against the single column the formula supplies, so
  # the formula fit really is over-identified and the controls are live.
  expect_identical(dim(psi(init)), c(2L, 200L))
  expect_grid_matches(
    overid_argument_sets,
    one_step,
    two_step,
    one_step(list(overid_maxiter = 200L))
  )
  expect_identical(
    one_step(list(overid_maxiter = 0L))@weight_matrix,
    diag(2)
  )
})

# ---- Phase 5: Broom tidiers -------------------------------------------------

test_that("tidy() returns expected columns", {
  m <- make_fitted_regression()
  td <- generics::tidy(m)
  expect_true(is.data.frame(td))
  expect_true(all(
    c("term", "estimate", "std.error", "statistic", "p.value", "s.value") %in%
      names(td)
  ))
  expect_equal(nrow(td), 3)
  expect_equal(td$term, c("intercept", "wt", "hp"))
})

test_that("tidy() with conf.int adds CI columns", {
  m <- make_fitted_regression()
  td <- generics::tidy(m, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in% names(td)))
})

test_that("tidy() conf.level works", {
  m <- make_fitted_regression()
  td_95 <- generics::tidy(m, conf.int = TRUE, conf.level = 0.95)
  td_99 <- generics::tidy(m, conf.int = TRUE, conf.level = 0.99)
  # 99% CI should be wider
  expect_true(
    td_99$conf.high[1] - td_99$conf.low[1] >
      td_95$conf.high[1] - td_95$conf.low[1]
  )
})

test_that("glance() returns expected columns", {
  m <- make_fitted_regression()
  gl <- generics::glance(m)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1)
  expect_true(all(c("nobs", "npar", "estimator") %in% names(gl)))
  expect_equal(gl$nobs, 32L)
  expect_equal(gl$npar, 3L)
  expect_equal(gl$estimator, "MEstimator")
})

test_that("tidy() works for GMMEstimator", {
  y <- c(1, 2, 3, 4, 5)
  g <- GMMEstimator(
    stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
    init = c(mean = 0)
  )
  g <- estimate(g)
  td <- generics::tidy(g)
  expect_equal(td$term, "mean")
  gl <- generics::glance(g)
  expect_equal(gl$estimator, "GMMEstimator")
})

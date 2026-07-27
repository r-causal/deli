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
# result relative to the fit without it. Once the subset and
# finite_correction formals are added to the wrappers this holds; until then
# the wrappers route
# these names into ... and on to the estimating equation, which rejects them
# with an unused-argument error, so the value assertions below error.

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
# m_estimate() and gmm_estimate() construct an estimator and call estimate() in
# a single call. Every constructor argument and every estimate() argument is
# handed straight through, so a one-step fit must be the same fit as the
# two-step form built from the same estimating function. The documentation
# depends on that equivalence: examples and articles use the one-step form
# everywhere the two-step form is not itself the subject. These tests hold the
# pass-through in place, and compare with expect_identical() rather than a
# tolerance so that a refactor which quietly reroutes an argument is caught
# instead of being absorbed by slack.

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

expect_same_fit <- function(one_step, two_step, label) {
  expect_identical(coef(one_step), coef(two_step), info = label)
  expect_identical(vcov(one_step), vcov(two_step), info = label)
}

# Linear regression on mtcars: three parameters, solved reliably from zero
# starting values by every solver under test.
mtcars_regression_psi <- function() {
  X <- stats::model.matrix(mpg ~ wt + hp, data = mtcars)
  y <- mtcars$mpg
  function(theta) ee_regression(theta, X = X, y = y, model = "linear")
}

mtcars_regression_init <- c(`(Intercept)` = 0, wt = 0, hp = 0)

# The argument surface the documentation actually passes, one set per argument
# plus a set carrying all of them at once.
m_argument_sets <- list(
  "defaults" = list(),
  "solver = NULL" = list(solver = NULL),
  "solver = rootSolve" = list(solver = "rootSolve"),
  "solver = lm" = list(solver = "lm"),
  "deriv_method = capprox" = list(deriv_method = "capprox"),
  "deriv_method = exact" = list(deriv_method = "exact"),
  "maxiter and tolerance" = list(maxiter = 200L, tolerance = 1e-7),
  "dx" = list(dx = 1e-6),
  "allow_pinv = FALSE" = list(allow_pinv = FALSE),
  "subset" = list(subset = 1L),
  "finite_correction = HC1" = list(finite_correction = "HC1"),
  "every argument at once" = list(
    solver = "lm",
    maxiter = 400L,
    tolerance = 1e-8,
    deriv_method = "exact",
    dx = 1e-7,
    allow_pinv = FALSE,
    subset = c(1L, 2L),
    finite_correction = "HC1"
  )
)

# GMM minimizes rather than root-finds, so the solver names are stats::optim
# methods, and it carries the two over-identification controls as well.
gmm_argument_sets <- list(
  "defaults" = list(),
  "solver = NULL" = list(solver = NULL),
  "solver = BFGS" = list(solver = "BFGS"),
  "deriv_method = capprox" = list(deriv_method = "capprox"),
  "deriv_method = exact" = list(deriv_method = "exact"),
  "maxiter and tolerance" = list(maxiter = 200L, tolerance = 1e-7),
  "dx" = list(dx = 1e-6),
  "allow_pinv = FALSE" = list(allow_pinv = FALSE),
  "subset" = list(subset = 1L),
  "finite_correction = HC1" = list(finite_correction = "HC1"),
  "overid controls" = list(overid_maxiter = 25L, overid_tolerance = 1e-7),
  "every argument at once" = list(
    solver = "BFGS",
    maxiter = 400L,
    tolerance = 1e-8,
    deriv_method = "exact",
    dx = 1e-7,
    allow_pinv = FALSE,
    subset = c(1L, 2L),
    finite_correction = "HC1",
    overid_maxiter = 25L,
    overid_tolerance = 1e-7
  )
)

test_that("m_estimate() function interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()

  for (label in names(m_argument_sets)) {
    args <- m_argument_sets[[label]]
    one_step <- do.call(
      m_estimate,
      c(
        list(stacked_equations = psi, init = mtcars_regression_init),
        args
      )
    )
    two_step <- fit_two_step(
      MEstimator,
      psi,
      mtcars_regression_init,
      args
    )
    expect_same_fit(one_step, two_step, label)
  }
})

test_that("m_estimate() formula interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()

  for (label in names(m_argument_sets)) {
    args <- m_argument_sets[[label]]
    one_step <- do.call(
      m_estimate,
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
    two_step <- fit_two_step(
      MEstimator,
      psi,
      mtcars_regression_init,
      args
    )
    expect_same_fit(one_step, two_step, label)
  }
})

test_that("gmm_estimate() function interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()

  for (label in names(gmm_argument_sets)) {
    args <- gmm_argument_sets[[label]]
    one_step <- do.call(
      gmm_estimate,
      c(
        list(stacked_equations = psi, init = mtcars_regression_init),
        args
      )
    )
    two_step <- fit_two_step(
      GMMEstimator,
      psi,
      mtcars_regression_init,
      args
    )
    expect_same_fit(one_step, two_step, label)
  }
})

test_that("gmm_estimate() formula interface reproduces the two-step fit", {
  psi <- mtcars_regression_psi()

  for (label in names(gmm_argument_sets)) {
    args <- gmm_argument_sets[[label]]
    one_step <- do.call(
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
    two_step <- fit_two_step(
      GMMEstimator,
      psi,
      mtcars_regression_init,
      args
    )
    expect_same_fit(one_step, two_step, label)
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

make_iv_psi <- function(d) {
  function(theta) {
    rbind(
      d$z1 * (d$y - theta[1] * d$a),
      d$z2 * (d$y - theta[1] * d$a)
    )
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

test_that("gmm_estimate() over-identified fit reproduces the two-step fit", {
  d <- make_iv_data()
  psi <- make_iv_psi(d)
  init <- c(effect = 0)

  overid_argument_sets <- list(
    "converged updating" = list(overid_maxiter = 200L),
    "no updating step" = list(overid_maxiter = 0L),
    "loose overid_tolerance" = list(
      overid_maxiter = 200L,
      overid_tolerance = 1
    )
  )

  fits <- list()
  for (label in names(overid_argument_sets)) {
    args <- overid_argument_sets[[label]]
    one_step <- do.call(
      gmm_estimate,
      c(list(stacked_equations = psi, init = init), args)
    )
    two_step <- fit_two_step(GMMEstimator, psi, init, args)
    expect_same_fit(one_step, two_step, label)
    fits[[label]] <- one_step
  }

  # Both controls change the answer, so the comparisons above are not vacuous.
  # Zero updating iterations leaves the identity weight matrix in place, and a
  # tolerance of 1 stops the updating after a single step.
  expect_identical(fits[["no updating step"]]@weight_matrix, diag(2))
  expect_false(identical(
    coef(fits[["converged updating"]]),
    coef(fits[["no updating step"]])
  ))
  expect_false(identical(
    coef(fits[["converged updating"]]),
    coef(fits[["loose overid_tolerance"]])
  ))
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

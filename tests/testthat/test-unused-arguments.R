# The entry points covered here are the ones that carry settings a wrong name
# could displace: `m_estimate()` and `gmm_estimate()` in their function form,
# `estimate()`, `delta_method()`, `confidence_bands()`, and the four inference
# generics. Each must reject an unrecognized name instead of discarding it. A
# swallowed name is worse than a hard error: `tolerence` silently keeps the
# default tolerance, and `deriv_methods` silently keeps the default Jacobian, so
# a typo changes a reported number without changing any visible code.
#
# The base-generic methods also carry a `...` they never forward, and they
# divide on whether a wrong name in it could displace anything. `summary()`,
# `confint()`, `generics::tidy()`, `print()`, and `model.matrix()` each carry a
# setting a typo could quietly redirect, so each of them is guarded the same way
# and the section after the inference generics covers them. `coef()`, `vcov()`,
# `nobs()`, and `generics::glance()` have no optional argument at all, so a
# stray name cannot change what they return and they stay tolerant, which is the
# convention for those generics.
#
# Each guarded entry point gets two tests. The first passes a name a user might
# plausibly reach for by mistake, either a misspelling or a name borrowed from a
# similar function, and requires the `rlib_error_dots_nonempty` condition that
# `rlang::check_dots_empty()` raises. The second is the positive control: every
# named argument of that interface at once, which fails if a guard is placed
# where a real argument would reach it.
#
# No wrong name used below is a prefix of any formal, because R matches a
# supplied name that is a prefix of exactly one formal to that formal. A dropped
# final letter (`n_draw` for `n_draws`) therefore never reaches `...` in the
# first place and cannot demonstrate anything. Pluralizing, shifting case, or
# substituting an interior letter leaves no such escape.
#
# The formula methods are deliberately excluded. Their `...` reaches the
# estimating equation, which rejects names it does not take, and the last
# section locks that forwarding down so a guard is never added there by
# analogy.

# ---- Helpers -----------------------------------------------------------------

# Mean and variance of a fixed sample: two parameters, so `subset` and the
# matrix-valued arguments have something to act on, and small enough that every
# solver reaches the root from the zero-ish start.
dots_psi <- function() {
  y <- c(1, 2, 3, 4, 5, 6, 7, 8)
  function(theta) rbind(y - theta[1], (y - theta[1])^2 - theta[2])
}

dots_init <- c(mean = 0, variance = 1)

dots_fit <- function() {
  m_estimate(stacked_equations = dots_psi(), init = dots_init)
}

# A custom solver satisfying the documented contract: takes `stacked_equations`
# and `init`, returns the solved vector. Minimizing the squared summed score
# finds the same root the default solver does.
score_solver <- function(stacked_equations, init) {
  stats::nlm(function(theta) sum(stacked_equations(theta)^2), p = init)$estimate
}

# The GMM path hands its solver a scalar objective rather than a score vector,
# so the custom solver minimizes it directly.
objective_solver <- function(stacked_equations, init) {
  stats::optim(par = init, fn = stacked_equations, method = "BFGS")$par
}

exp_first <- function(theta) exp(theta[1])

# Evaluates a quoted call that must trip the dots guard, and asserts that the
# error is attributed to that same call. The call is quoted rather than passed
# as an expression so the expectation can compare against it directly, and
# `eval_bare()` runs it in the caller's environment without interposing a frame
# of its own that could change which frame the guard reports.
expect_reported_call <- function(call, env = rlang::caller_env()) {
  cnd <- rlang::catch_cnd(
    rlang::eval_bare(call, env),
    classes = "rlib_error_dots_nonempty"
  )
  expect_equal(conditionCall(cnd), call)
}

# ---- m_estimate() ------------------------------------------------------------

test_that("m_estimate() rejects a misspelled argument in the function form", {
  expect_error(
    m_estimate(
      stacked_equations = dots_psi(),
      init = dots_init,
      tolerence = 1e-12
    ),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("m_estimate() rejects an unnamed second argument", {
  # `init` follows `...` in the signature, so an unnamed second argument lands
  # in `...` rather than in `init`. The guard names the offending position and
  # asks whether a name was forgotten, which is the actionable form of what was
  # previously reported as a missing `init`.
  expect_error(
    m_estimate(dots_psi(), dots_init),
    "Did you forget to name an argument",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("m_estimate() accepts every argument of the function interface", {
  fit <- m_estimate(
    stacked_equations = dots_psi(),
    init = dots_init,
    subset = c(1L, 2L),
    finite_correction = "HC1",
    solver = "rootSolve",
    maxiter = 500L,
    tolerance = 1e-10,
    deriv_method = "fapprox",
    dx = 1e-7,
    allow_pinv = TRUE
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-6)
})

test_that("m_estimate() accepts a function passed as solver", {
  fit <- m_estimate(
    stacked_equations = dots_psi(),
    init = dots_init,
    solver = score_solver
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-6)
})

# ---- gmm_estimate() ----------------------------------------------------------

test_that("gmm_estimate() rejects a misspelled argument in the function form", {
  expect_error(
    gmm_estimate(
      stacked_equations = dots_psi(),
      init = dots_init,
      maxIter = 100L
    ),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("gmm_estimate() rejects an unnamed second argument", {
  expect_error(
    gmm_estimate(dots_psi(), dots_init),
    "Did you forget to name an argument",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("gmm_estimate() accepts every argument of the function interface", {
  fit <- gmm_estimate(
    stacked_equations = dots_psi(),
    init = dots_init,
    subset = c(1L, 2L),
    finite_correction = "HC1",
    solver = "BFGS",
    maxiter = 5000L,
    tolerance = 1e-10,
    deriv_method = "fapprox",
    dx = 1e-7,
    allow_pinv = TRUE,
    overid_maxiter = 10L,
    overid_tolerance = 1e-9
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-3)
})

test_that("gmm_estimate() accepts a function passed as solver", {
  fit <- gmm_estimate(
    stacked_equations = dots_psi(),
    init = dots_init,
    solver = objective_solver
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-3)
})

# ---- estimate() --------------------------------------------------------------

test_that("estimate() rejects a misspelled argument for an MEstimator", {
  obj <- MEstimator(stacked_equations = dots_psi(), init = dots_init)
  expect_error(
    estimate(obj, tolerence = 1e-12),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("estimate() rejects a misspelled argument for a GMMEstimator", {
  obj <- GMMEstimator(stacked_equations = dots_psi(), init = dots_init)
  expect_error(
    estimate(obj, maxIter = 100L),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("estimate() accepts every documented argument for an MEstimator", {
  obj <- MEstimator(stacked_equations = dots_psi(), init = dots_init)
  fit <- estimate(
    obj,
    solver = "rootSolve",
    maxiter = 500L,
    tolerance = 1e-10,
    deriv_method = "fapprox",
    dx = 1e-7,
    allow_pinv = TRUE
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-6)
})

test_that("estimate() accepts every documented argument for a GMMEstimator", {
  obj <- GMMEstimator(stacked_equations = dots_psi(), init = dots_init)
  fit <- estimate(
    obj,
    solver = "BFGS",
    maxiter = 5000L,
    tolerance = 1e-10,
    deriv_method = "fapprox",
    dx = 1e-7,
    allow_pinv = TRUE
  )
  expect_equal(unname(coef(fit)[1]), 4.5, tolerance = 1e-3)
})

test_that("estimate() accepts a function passed as solver for both classes", {
  m <- estimate(
    MEstimator(stacked_equations = dots_psi(), init = dots_init),
    solver = score_solver
  )
  g <- estimate(
    GMMEstimator(stacked_equations = dots_psi(), init = dots_init),
    solver = objective_solver
  )
  expect_equal(unname(coef(m)[1]), 4.5, tolerance = 1e-6)
  expect_equal(unname(coef(g)[1]), 4.5, tolerance = 1e-3)
})

# ---- delta_method() ----------------------------------------------------------
#
# The worst case of the defect class: `deriv_methods` does not change a setting
# that a later check would catch, it quietly substitutes a finite-difference
# Jacobian for the exact one and changes the returned variance.

test_that("delta_method() rejects a misspelled argument for an estimator", {
  fit <- dots_fit()
  expect_error(
    delta_method(fit, transform = exp_first, deriv_methods = "exact"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("delta_method() rejects a misspelled argument for a numeric vector", {
  fit <- dots_fit()
  expect_error(
    delta_method(
      coef(fit),
      transform = exp_first,
      covariance = vcov(fit),
      deriv_methods = "exact"
    ),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("delta_method() accepts every documented argument for an estimator", {
  fit <- dots_fit()
  cov_exact <- delta_method(
    fit,
    transform = exp_first,
    covariance = NULL,
    deriv_method = "exact",
    dx = 1e-7
  )
  expect_equal(dim(cov_exact), c(1L, 1L))
  expect_true(cov_exact[1, 1] > 0)
})

test_that("delta_method() accepts every documented argument for a vector", {
  fit <- dots_fit()
  cov_exact <- delta_method(
    coef(fit),
    transform = exp_first,
    covariance = vcov(fit),
    deriv_method = "exact",
    dx = 1e-7
  )
  expect_equal(dim(cov_exact), c(1L, 1L))
  expect_true(cov_exact[1, 1] > 0)
})

# ---- confidence_bands() ------------------------------------------------------
#
# The Bonferroni method keeps these tests free of the multivariate-normal draw,
# so `n_draws` and `seed` are checked as accepted formals rather than for their
# effect, which test-confidence-bands.R already covers.

test_that("confidence_bands() rejects a misspelled argument for an estimator", {
  fit <- dots_fit()
  expect_error(
    confidence_bands(fit, methods = "bonferroni"),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("confidence_bands() rejects a misspelled argument for a vector", {
  fit <- dots_fit()
  expect_error(
    confidence_bands(
      coef(fit),
      covariance = vcov(fit),
      methods = "bonferroni"
    ),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("confidence_bands() accepts every documented argument", {
  fit <- dots_fit()
  from_fit <- confidence_bands(
    fit,
    alpha = 0.1,
    method = "bonferroni",
    n_draws = 1000L,
    seed = 1,
    subset = 1L,
    covariance = NULL
  )
  from_vector <- confidence_bands(
    coef(fit),
    alpha = 0.1,
    method = "bonferroni",
    n_draws = 1000L,
    seed = 1,
    subset = 1L,
    covariance = vcov(fit)
  )
  expect_equal(dim(from_fit), c(1L, 2L))
  expect_equal(unname(from_vector), unname(from_fit))
})

# ---- Inference generics ------------------------------------------------------

test_that("confidence_intervals() rejects a misspelled argument", {
  fit <- dots_fit()
  expect_error(
    confidence_intervals(fit, level = 0.95),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("z_scores(), p_values(), and s_values() reject a misspelled name", {
  fit <- dots_fit()
  expect_error(z_scores(fit, nulls = 1), class = "rlib_error_dots_nonempty")
  expect_error(p_values(fit, nulls = 1), class = "rlib_error_dots_nonempty")
  expect_error(s_values(fit, nulls = 1), class = "rlib_error_dots_nonempty")
})

test_that("influence_functions() rejects a misspelled argument", {
  fit <- dots_fit()
  expect_error(
    influence_functions(fit, allow_pinverse = TRUE),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("the inference generics accept every documented argument", {
  fit <- dots_fit()
  expect_equal(dim(confidence_intervals(fit, alpha = 0.1)), c(2L, 2L))
  expect_length(z_scores(fit, null = c(4, 5)), 2L)
  expect_length(p_values(fit, null = c(4, 5)), 2L)
  expect_length(s_values(fit, null = c(4, 5)), 2L)
  expect_equal(dim(influence_functions(fit, allow_pinv = FALSE)), c(8L, 2L))
})

# ---- Base-generic methods that carry a displaceable setting ------------------
#
# `summary()`, `confint()`, `generics::tidy()`, `print()`, and `model.matrix()`
# each have an optional argument a wrong name can displace: the significance
# level and the displayed subset, the interval level and the parameters chosen,
# whether an interval is returned and at what level, the displayed subset again,
# and the data a design is built from. Each wrong name changes a reported number
# or a reported row set while leaving the call looking like the one that was
# meant, which is the failure the first section of this file describes, so each
# of the five carries the same guard.
#
# `coef()`, `vcov()`, `nobs()`, and `generics::glance()` have no optional
# argument for a name to displace, so a stray name cannot change what they
# return. They stay tolerant, and the last test here pins that.
#
# Which wrong name reaches `...` depends on where the argument it was meant for
# sits in the signature. `alpha` and `subset` precede `summary()`'s dots, so
# `alph`, `a`, `subse`, and `s` all resolve by partial matching before the dots
# collect anything and no guard can see them. The same holds for `lev` and `par`
# on `confint()`, `conf.i` and `conf.l` on `tidy()`, and `dat` on
# `model.matrix()`; a bare `conf` on `tidy()` matches two formals and R rejects
# it on its own. `print()` is the exception, because `subset` follows its dots
# and therefore matches only exactly: every misspelling of it, `subse` and `sub`
# included, reaches the guard.

dots_gmm_fit <- function() {
  gmm_estimate(stacked_equations = dots_psi(), init = dots_init)
}

# `model.matrix()` reports on a recorded design, which only the formula
# interface records, so its guard is checked on a formula fit. The factor is
# what makes a `contrasts.arg` a request the method could conceivably have
# honored rather than a name with nothing to act on.
dots_design_data <- function() {
  data <- mtcars
  data$gear <- factor(data$gear)
  data
}

dots_design_fit <- function(data = dots_design_data()) {
  m_estimate(
    mpg ~ wt + gear,
    data = data,
    .ee = ee_regression,
    model = "linear"
  )
}

test_that("summary() rejects a misspelled argument", {
  fit <- dots_fit()
  expect_error(
    summary(fit, alfa = 0.1),
    "alfa",
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    summary(fit, subst = 1),
    "subst",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("summary() rejects a misspelled argument for a GMMEstimator", {
  fit <- dots_gmm_fit()
  expect_error(
    summary(fit, alfa = 0.1),
    "alfa",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("summary() accepts both of its own arguments", {
  fit <- dots_fit()
  s <- summary(fit, alpha = 0.1, subset = 1L)
  expect_equal(s@alpha, 0.1)
  expect_named(s@theta, "mean")
})

test_that("confint() rejects a misspelled level", {
  fit <- dots_fit()
  expect_error(
    confint(fit, levl = 0.9),
    "levl",
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    confint(fit, lvl = 0.9),
    "lvl",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("confint() rejects a misspelled level for a GMMEstimator", {
  fit <- dots_gmm_fit()
  expect_error(
    confint(fit, levl = 0.9),
    "levl",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("confint() accepts both of its own arguments", {
  fit <- dots_fit()
  wide <- confint(fit, level = 0.99)
  narrow <- confint(fit, level = 0.5)
  expect_true(all(
    wide[, "upper"] - wide[, "lower"] > narrow[, "upper"] - narrow[, "lower"]
  ))
  expect_identical(rownames(confint(fit, parm = "variance")), "variance")
  expect_equal(
    confint(fit, parm = 1L, level = 0.5),
    narrow[1L, , drop = FALSE]
  )
})

test_that("tidy() rejects a misspelled argument", {
  fit <- dots_fit()
  expect_error(
    generics::tidy(fit, conf.intt = TRUE),
    "conf.intt",
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    generics::tidy(fit, conf.int = TRUE, conf.levell = 0.5),
    "conf.levell",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("tidy() accepts both of its own arguments", {
  fit <- dots_fit()
  plain <- generics::tidy(fit)
  narrow <- generics::tidy(fit, conf.int = TRUE, conf.level = 0.5)
  expect_false("conf.low" %in% names(plain))
  expect_true(all(c("conf.low", "conf.high") %in% names(narrow)))
  expect_equal(
    narrow$conf.low,
    unname(confint(fit, level = 0.5)[, "lower"])
  )
})

test_that("print() rejects a misspelled subset", {
  fit <- dots_fit()
  expect_error(
    print(fit, subst = 1),
    "subst",
    class = "rlib_error_dots_nonempty"
  )
  # `subset` follows the dots, so it matches only exactly and a dropped final
  # letter reaches the guard rather than the argument it was meant for.
  expect_error(
    print(fit, subse = 1),
    "subse",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("print() accepts its own subset", {
  fit <- dots_fit()
  shown <- testthat::capture_messages(print(fit, subset = 1L))
  expect_true(any(grepl("mean", shown, fixed = TRUE)))
  expect_false(any(grepl("variance", shown, fixed = TRUE)))
})

test_that("model.matrix() rejects an argument it has no answer for", {
  fit <- dots_design_fit()
  # A design is coded with the contrasts the fit recorded rather than with a
  # coding supplied at the call, so there is no `contrasts.arg` to honor and a
  # supplied one has to be refused rather than dropped.
  expect_error(
    stats::model.matrix(fit, contrasts.arg = list(gear = "contr.sum")),
    "contrasts.arg",
    class = "rlib_error_dots_nonempty"
  )
})

test_that("model.matrix() accepts its own data", {
  data <- dots_design_data()
  fit <- dots_design_fit(data)
  expect_identical(stats::model.matrix(fit), fit@model_spec$X)
  expect_identical(
    nrow(stats::model.matrix(fit, data = data[1:3, ])),
    3L
  )
})

test_that("coef(), vcov(), nobs(), and glance() stay tolerant of a stray name", {
  fit <- dots_fit()
  expect_equal(unname(coef(fit, junk = 1)[1]), 4.5, tolerance = 1e-6)
  expect_equal(dim(vcov(fit, junk = 1)), c(2L, 2L))
  expect_equal(nobs(fit, junk = 1), 8L)
  expect_equal(generics::glance(fit, junk = 1)$nobs, 8L)
})

test_that("the dots guard names the user's call at each guarded method", {
  fit <- dots_fit()
  design_fit <- dots_design_fit()

  expect_reported_call(quote(summary(fit, alfa = 0.1)))
  expect_reported_call(quote(confint(fit, levl = 0.9)))
  expect_reported_call(quote(generics::tidy(fit, conf.intt = TRUE)))
  expect_reported_call(quote(print(fit, subst = 1)))
  expect_reported_call(quote(stats::model.matrix(
    design_fit,
    contrasts.arg = list(gear = "contr.sum")
  )))
})

# ---- Guard ordering ----------------------------------------------------------
#
# Five methods run both `rlang::check_dots_empty()` and `check_estimated()`:
# `confidence_intervals()`, `z_scores()`, `influence_functions()`, and the
# estimator methods of `confidence_bands()` and `delta_method()`. Each runs the
# dots guard first. A wrong name is a mistake the user can see in the call they
# just wrote, so it is the more useful of the two errors to report. Reordering
# the two checks would silently change which one a user gets. `p_values()` and
# `s_values()` reach `check_estimated()` only through `z_scores()`, so their
# ordering follows from the one asserted here.

test_that("the dots guard reports before the unfitted-object check", {
  obj <- MEstimator(stacked_equations = dots_psi(), init = dots_init)

  expect_error(
    confidence_intervals(obj, level = 0.95),
    class = "rlib_error_dots_nonempty"
  )
  expect_error(z_scores(obj, nulls = 1), class = "rlib_error_dots_nonempty")
  expect_error(
    influence_functions(obj, allow_pinverse = TRUE),
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    confidence_bands(obj, methods = "bonferroni"),
    class = "rlib_error_dots_nonempty"
  )
  expect_error(
    delta_method(obj, transform = exp_first, deriv_methods = "exact"),
    class = "rlib_error_dots_nonempty"
  )

  # Without the wrong name the same unfitted object reaches check_estimated(),
  # so the assertions above are about ordering and not about the object.
  unfitted <- "Cannot compute inference before"
  expect_error(confidence_intervals(obj), unfitted)
  expect_error(z_scores(obj), unfitted)
  expect_error(influence_functions(obj), unfitted)
  expect_error(confidence_bands(obj), unfitted)
  expect_error(delta_method(obj, transform = exp_first), unfitted)
})

# ---- The call the guard reports ----------------------------------------------
#
# A guard whose only job is to point at a typo has to point at the call that
# contains it. Inside an S7 method the default `call = caller_env()` of
# `rlang::check_dots_empty()` resolves to the dispatched method, so the error
# header reads `method(confidence_bands, class_double)`, which is not a callable
# expression and appears nowhere in the user's code. Passing
# `call = rlang::caller_env()` from the method body resolves to the generic's
# frame instead, whose call is the one the user wrote. The two `UseMethod()`
# wrappers need no argument: their caller is already `m_estimate()` or
# `gmm_estimate()`.

test_that("the dots guard names the user's call at every S7 entry point", {
  fit <- dots_fit()
  m_obj <- MEstimator(stacked_equations = dots_psi(), init = dots_init)
  g_obj <- GMMEstimator(stacked_equations = dots_psi(), init = dots_init)
  theta <- coef(fit)
  covariance <- vcov(fit)

  expect_reported_call(quote(estimate(m_obj, tolerence = 1e-12)))
  expect_reported_call(quote(estimate(g_obj, maxIter = 100L)))
  expect_reported_call(quote(confidence_intervals(fit, level = 0.95)))
  expect_reported_call(quote(z_scores(fit, nulls = 1)))
  expect_reported_call(quote(p_values(fit, nulls = 1)))
  expect_reported_call(quote(s_values(fit, nulls = 1)))
  expect_reported_call(quote(influence_functions(fit, allow_pinverse = TRUE)))
  expect_reported_call(quote(confidence_bands(fit, methods = "bonferroni")))
  expect_reported_call(quote(confidence_bands(
    theta,
    covariance = covariance,
    methods = "bonferroni"
  )))
  expect_reported_call(quote(delta_method(
    fit,
    transform = exp_first,
    deriv_methods = "exact"
  )))
  expect_reported_call(quote(delta_method(
    theta,
    transform = exp_first,
    covariance = covariance,
    deriv_methods = "exact"
  )))
})

test_that("the dots guard names the user's call at both S3 wrappers", {
  psi <- dots_psi()

  expect_reported_call(quote(m_estimate(
    stacked_equations = psi,
    init = dots_init,
    tolerence = 1e-12
  )))
  expect_reported_call(quote(gmm_estimate(
    stacked_equations = psi,
    init = dots_init,
    maxIter = 100L
  )))
})

# ---- The formula methods keep forwarding dots --------------------------------
#
# `m_estimate.formula()` and `gmm_estimate.formula()` route `...` to the
# estimating equation, evaluated against `data`. That is the documented
# contract, and it is the reason no dots guard belongs in those two methods.

test_that("m_estimate() formula dots reach the estimating equation", {
  set.seed(4)
  d <- mtcars
  d$w <- stats::runif(nrow(d), 0.5, 1.5)

  weighted <- m_estimate(
    mpg ~ wt + hp,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
  unweighted <- m_estimate(
    mpg ~ wt + hp,
    data = d,
    .ee = ee_regression,
    model = "linear"
  )
  oracle <- stats::lm(mpg ~ wt + hp, data = d, weights = w)

  expect_equal(unname(coef(weighted)), unname(coef(oracle)), tolerance = 1e-6)
  # Without this the equality above would hold whether or not `weights` arrived.
  expect_false(isTRUE(all.equal(coef(weighted), coef(unweighted))))
})

test_that("gmm_estimate() formula dots reach the estimating equation", {
  set.seed(4)
  d <- mtcars
  d$w <- stats::runif(nrow(d), 0.5, 1.5)

  weighted <- gmm_estimate(
    mpg ~ wt + hp,
    data = d,
    .ee = ee_regression,
    model = "linear",
    weights = w
  )
  oracle <- stats::lm(mpg ~ wt + hp, data = d, weights = w)

  expect_equal(unname(coef(weighted)), unname(coef(oracle)), tolerance = 1e-5)
})

test_that("the formula interface rejects a name the .ee does not take", {
  expect_error(
    m_estimate(
      mpg ~ wt + hp,
      data = mtcars,
      .ee = ee_regression,
      model = "linear",
      init = c(`(Intercept)` = 0, wt = 0, hp = 0),
      tolerence = 1e-12
    ),
    "unused argument"
  )
})

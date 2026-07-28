# Each deli utility function mirrors a base R primitive, and the reference
# documentation names the counterpart and says where each may be used. The
# load-bearing half of that claim is the exact-mode split. `deriv_method =
# "exact"` is forward-mode automatic differentiation dispatched through the `Ops`
# and `Math` group generics, so a compiled base R primitive outside both groups
# cannot accept a tangent-carrying argument, while the deli version can. Five of
# the seven pairs split that way. `digamma` is a `Math` group member with a
# tangent rule and `identity` is pure pass-through, so those two work under every
# derivative method.
#
#   deli function            base R counterpart              base R under exact
#   ------------------------ ------------------------------- ------------------
#   inverse_logit(x)         stats::plogis(x)                 errors
#   logit(p)                 stats::qlogis(p)                 errors
#   standard_normal_cdf(x)   stats::pnorm(x)                  errors
#   standard_normal_pdf(x)   stats::dnorm(x)                  errors
#   deli_polygamma(n, x)     base::psigamma(x, deriv = n)     errors
#   deli_digamma(z)          base::digamma(z)                 works
#   identity_transform(v)    base::identity(v)                works
#
# This file is the regression lock for that split, in both of the positions the
# documented rule covers: inside a transform passed to `delta_method()` and
# inside an estimating equation. If the autodiff layer later learns to handle
# `plogis`, or the `Math` group rules change, the documentation becomes wrong and
# a failure here is what says so.
#
# The base R primitives raise "Non-numeric argument to mathematical function"
# from compiled code, with capitalization that varies between primitives. deli
# catches that error inside auto_differentiation(), reads the offending function
# off the failing call, and re-raises its own abort naming both the function and
# its deli replacement, with the class deli_exact_unsupported_function. These
# tests therefore assert on the class and on the two function names, which is
# what the documentation promises, rather than on the surrounding prose.
#
# Numeric agreement between each deli function and its counterpart on plain
# numeric input lives with the functions themselves, in test-transforms.R and
# test-math-utils.R, alongside the pole and argument-order assertions for
# deli_digamma() and deli_polygamma().

# A fitted logistic regression. Its coefficients are on the log-odds scale, which
# the scale helpers below map to whatever domain a given pair needs.
mtcars_logistic_fit <- function() {
  m_estimate(
    vs ~ mpg,
    data = mtcars,
    .ee = ee_regression,
    model = "logistic"
  )
}

# Scale helpers shared by both halves of every pair. Each is plain arithmetic or
# a `Math` group member, so it carries a tangent under exact mode whichever outer
# function is applied. The deli transform and its base R partner therefore differ
# in exactly one place, the outer call, which is what makes a failure
# attributable to the function under test.
log_odds_at_20 <- function(theta) theta[1] + theta[2] * 20
probability_at_20 <- function(theta) 1 / (1 + exp(-log_odds_at_20(theta)))
positive_scale <- function(theta) exp(theta[2])

# Assert the documented behavior for one row of the correspondence table:
#
#   * the deli transform differentiates under exact mode, giving a strictly
#     positive variance that reproduces the central-difference reference;
#   * the base R transform is fine under `capprox` and agrees numerically there,
#     so exact mode is the discriminator rather than the function itself;
#   * under exact mode the base R transform either errors or matches the deli
#     result, according to `base_exact_works`.
#
# When it errors, the abort must carry the deli_exact_unsupported_function class
# and name both halves of the pair, so a user reading it learns which call
# failed and what to write instead. Both names come from `pair`, which is
# written as "<deli function> / <base R function>".
#
# The strictly positive variance keeps the comparisons from passing vacuously
# against a zero Jacobian. The 1e-5 tolerance follows the package convention for
# comparing an exact Jacobian against a central difference, which carries roughly
# 1e-6 numerical error at the default 1e-9 step size.
expect_exact_mode_split <- function(
  deli_transform,
  base_transform,
  pair,
  base_exact_works = FALSE
) {
  fit <- mtcars_logistic_fit()
  deli_exact <- delta_method(
    fit,
    transform = deli_transform,
    deriv_method = "exact"
  )
  deli_capprox <- delta_method(
    fit,
    transform = deli_transform,
    deriv_method = "capprox"
  )
  base_capprox <- delta_method(
    fit,
    transform = base_transform,
    deriv_method = "capprox"
  )

  testthat::expect_true(
    all(diag(deli_exact) > 0),
    label = paste0(pair, ": deli exact variance is strictly positive")
  )
  testthat::expect_equal(
    deli_exact,
    deli_capprox,
    tolerance = 1e-5,
    label = paste0(pair, ": deli exact vs deli capprox")
  )
  testthat::expect_equal(
    base_capprox,
    deli_capprox,
    tolerance = 1e-5,
    label = paste0(pair, ": base R capprox vs deli capprox")
  )

  if (base_exact_works) {
    testthat::expect_equal(
      delta_method(fit, transform = base_transform, deriv_method = "exact"),
      deli_exact,
      label = paste0(pair, ": base R exact vs deli exact")
    )
  } else {
    err <- testthat::expect_error(
      delta_method(fit, transform = base_transform, deriv_method = "exact"),
      class = "deli_exact_unsupported_function"
    )
    # cli wraps the bullets, so runs of whitespace collapse to a single space
    # before the two names are matched.
    pair_names <- strsplit(pair, " / ", fixed = TRUE)[[1]]
    flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
    testthat::expect_match(flat, pair_names[2], fixed = TRUE)
    testthat::expect_match(flat, pair_names[1], fixed = TRUE)
  }
  invisible(deli_exact)
}

# ---- the five pairs where base R errors under exact mode --------------------

test_that("inverse_logit() differentiates under exact mode where plogis() errors", {
  expect_exact_mode_split(
    deli_transform = function(theta) inverse_logit(log_odds_at_20(theta)),
    base_transform = function(theta) plogis(log_odds_at_20(theta)),
    pair = "inverse_logit / plogis"
  )
})

test_that("logit() differentiates under exact mode where qlogis() errors", {
  expect_exact_mode_split(
    deli_transform = function(theta) logit(probability_at_20(theta)),
    base_transform = function(theta) qlogis(probability_at_20(theta)),
    pair = "logit / qlogis"
  )
})

test_that("standard_normal_cdf() differentiates under exact mode where pnorm() errors", {
  expect_exact_mode_split(
    deli_transform = function(theta) standard_normal_cdf(theta[2]),
    base_transform = function(theta) pnorm(theta[2]),
    pair = "standard_normal_cdf / pnorm"
  )
})

test_that("standard_normal_pdf() differentiates under exact mode where dnorm() errors", {
  expect_exact_mode_split(
    deli_transform = function(theta) standard_normal_pdf(theta[2]),
    base_transform = function(theta) dnorm(theta[2]),
    pair = "standard_normal_pdf / dnorm"
  )
})

test_that("deli_polygamma() differentiates under exact mode where psigamma() errors", {
  expect_exact_mode_split(
    deli_transform = function(theta) deli_polygamma(1, positive_scale(theta)),
    base_transform = function(theta) {
      psigamma(positive_scale(theta), deriv = 1)
    },
    pair = "deli_polygamma / psigamma"
  )
})

# ---- the two pairs where base R works under exact mode ----------------------

test_that("deli_digamma() and base::digamma() both differentiate under exact mode", {
  expect_exact_mode_split(
    deli_transform = function(theta) deli_digamma(positive_scale(theta)),
    base_transform = function(theta) digamma(positive_scale(theta)),
    pair = "deli_digamma / digamma",
    base_exact_works = TRUE
  )
})

test_that("identity_transform() and base::identity() both differentiate under exact mode", {
  expect_exact_mode_split(
    deli_transform = function(theta) identity_transform(theta[2]),
    base_transform = function(theta) identity(theta[2]),
    pair = "identity_transform / identity",
    base_exact_works = TRUE
  )
})

# ---- the same split inside an estimating equation ---------------------------
#
# The documented rule covers estimating equations as well as delta-method
# transforms, so the split is pinned in `psi` position too. An intercept-only
# score Y - g(theta) is the smallest equation that puts the transform where the
# bread has to differentiate it: root-finding runs on plain numeric theta, so a
# base R `psi` solves and then fails when the exact bread is built. The two
# functions covered here are the ones that appear in `psi` position in practice,
# a logistic mean and a probit mean; the delta-method tests above cover all seven
# pairs.

test_that("inverse_logit() in a psi solves under exact mode where plogis() errors", {
  set.seed(20240115)
  y <- rbinom(200, 1, 0.4)
  psi_deli <- function(theta) matrix(y - inverse_logit(theta[1]), nrow = 1)
  psi_base <- function(theta) matrix(y - plogis(theta[1]), nrow = 1)

  m <- m_estimate(
    stacked_equations = psi_deli,
    init = 0,
    deriv_method = "exact"
  )
  # The intercept-only logistic score solves at the log-odds of the sample mean.
  expect_equal(unname(coef(m)), logit(mean(y)))
  expect_gt(unname(sqrt(diag(vcov(m)))), 0)

  # The base R psi is fine under the finite-difference bread and gives the same
  # answer, so exact mode is the discriminator rather than the function itself.
  m_base <- m_estimate(
    stacked_equations = psi_base,
    init = 0,
    deriv_method = "capprox"
  )
  expect_equal(coef(m_base), coef(m))
  expect_equal(vcov(m_base), vcov(m), tolerance = 1e-6)

  expect_error(
    m_estimate(stacked_equations = psi_base, init = 0, deriv_method = "exact"),
    class = "deli_exact_unsupported_function"
  )
})

test_that("standard_normal_cdf() in a psi solves under exact mode where pnorm() errors", {
  set.seed(20240115)
  y <- rbinom(200, 1, 0.4)
  psi_deli <- function(theta) {
    matrix(y - standard_normal_cdf(theta[1]), nrow = 1)
  }
  psi_base <- function(theta) matrix(y - pnorm(theta[1]), nrow = 1)

  m <- m_estimate(
    stacked_equations = psi_deli,
    init = 0,
    deriv_method = "exact"
  )
  # The intercept-only probit score solves at the normal quantile of the sample
  # mean.
  expect_equal(unname(coef(m)), qnorm(mean(y)))
  expect_gt(unname(sqrt(diag(vcov(m)))), 0)

  m_base <- m_estimate(
    stacked_equations = psi_base,
    init = 0,
    deriv_method = "capprox"
  )
  expect_equal(coef(m_base), coef(m))
  expect_equal(vcov(m_base), vcov(m), tolerance = 1e-6)

  expect_error(
    m_estimate(stacked_equations = psi_base, init = 0, deriv_method = "exact"),
    class = "deli_exact_unsupported_function"
  )
})

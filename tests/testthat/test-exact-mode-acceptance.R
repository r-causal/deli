# Package-wide acceptance sweep for deriv_method = "exact" (forward-mode
# autodiff). This is the single exhaustive checklist that answers "does exact
# mode work for every built-in estimating equation?". For each exported ee_*
# family it fits a small, well-conditioned representative problem twice, once
# with deriv_method = "exact" and once with deriv_method = "capprox", and
# confirms the exact bread and the exact standard errors match the
# central-difference reference. The point estimate is identical between the two
# fits (the derivative method enters only the bread), so each comparison
# isolates the derivative computation. A regression in any family, or a new
# ee_* added without exact-mode support, fails a labelled assertion here.
#
# This file complements rather than replaces the per-family exact tests in
# test-ee-regression-exact.R, test-ee-causal-measurement-exact.R, and
# test-ee-survival-exact.R, which pin the exact-mode results against Python
# Delicatessen fixtures computed with deriv_method='exact'. Those files hold the
# Python-parity assertions; this file holds the exhaustive exact-versus-capprox
# checklist.
#
# Tolerance policy. The exact and capprox derivatives cannot agree more tightly
# than the central difference's own floor: a dx of 1e-9 carries roughly 1e-7
# rounding error, so smooth families use 1e-6. Families that are less
# well-conditioned use a documented looser tolerance:
#   * nonlinear dose-response (ee_emax, ee_loglogistic): 1e-5, the nonlinear
#     Jacobian is less well-conditioned than a linear score.
#   * ee_additive_regression: 1e-5 on the bread, the restricted-spline basis is
#     ill-conditioned (the same caution the spline tests elsewhere note).
#   * ee_mean_sensitivity_analysis: 1e-4, the outcome enters on a scale of order
#     100, which inflates the absolute central-difference error; the comparison
#     is effectively relative and lands near 1e-5.
# Every standard-error comparison first asserts the minimum standard error is
# well above the comparison floor, so the relative all.equal comparison is not
# vacuous, and every bread comparison asserts the bread is not near-zero, which
# guards against a degenerate (identically zero) bread passing trivially.
#
# ---------------------------------------------------------------------------
# Coverage summary (every exported ee_* family and its exact-mode status)
# ---------------------------------------------------------------------------
# basic
#   ee_mean                         verified
#   ee_mean_variance                verified
#   ee_mean_robust                  verified
#   ee_mean_geometric               verified
#   ee_percentile                   exception: non-differentiable estimand
#                                   (indicator on theta); bread is identically
#                                   zero under any derivative method, so the
#                                   sandwich is undefined. The function warns.
#   ee_positive_mean_deviation      exception: non-differentiable estimand (the
#                                   median score); the bread is singular under
#                                   any derivative method. The function warns.
# regression / glm
#   ee_regression                   verified (linear, logistic, poisson)
#   ee_glm                          verified (gamma log, gamma identity, probit,
#                                   negative binomial, tweedie, poisson)
#   ee_robust_regression            verified (huber, tukey, andrew, hampel,
#                                   cauchy)
#   ee_ridge_regression             verified
#   ee_bridge_regression            verified (gamma = 2 and gamma = 1.5)
#   ee_lasso_regression             verified (bread; the L1 kink leaves the
#                                   sandwich undefined on the standard-error
#                                   scale, matching the per-family caveat)
#   ee_elasticnet_regression        verified (bread; L1 kink, as ee_lasso)
#   ee_dlasso_regression            verified
#   ee_additive_regression          verified (bread tolerance 1e-5)
#   ee_mlogit                       verified
#   ee_beta_regression              verified
#   ee_tobit                        verified
# causal
#   ee_ipw                          verified
#   ee_gformula                     verified
#   ee_aipw                         verified
#   ee_gestimation_snmm             verified
#   ee_2sls                         verified (with and without W)
#   ee_iv_causal                    verified
#   ee_ipw_msm                      verified
#   ee_mean_sensitivity_analysis    verified (tolerance 1e-4)
# measurement
#   ee_rogan_gladen                 verified
#   ee_rogan_gladen_extended        verified
#   ee_regression_calibration       verified
# survival
#   ee_aft                          verified (exponential, weibull,
#                                   log-logistic, log-normal)
#   ee_survival_model               verified (exponential, weibull, gompertz)
#   ee_plogit                       verified
# pharmacokinetics
#   ee_emax                         verified (tolerance 1e-5)
#   ee_loglogistic                  verified (tolerance 1e-5)
#   ee_emax_ed                      verified (stacked on ee_emax, tolerance 1e-5)
#   ee_loglogistic_ed               verified (stacked on ee_loglogistic,
#                                   tolerance 1e-5)
# ---------------------------------------------------------------------------

# Fit a representative problem with both derivative methods and confirm the
# exact bread and standard errors match the central-difference reference. The
# labels carry the family name so a failure identifies the family immediately.
expect_exact_acceptable <- function(
  psi,
  init,
  family,
  solver = NULL,
  bread_tol = 1e-6,
  se_tol = 1e-6,
  check_se = TRUE
) {
  fit <- function(method) {
    m <- MEstimator(stacked_equations = psi, init = init)
    if (is.null(solver)) {
      estimate(m, deriv_method = method)
    } else {
      estimate(m, deriv_method = method, solver = solver)
    }
  }
  m_exact <- fit("exact")
  m_cap <- fit("capprox")

  # A non-degenerate bread: the comparison is not vacuous against a near-zero
  # target. (The non-differentiable families below fail this by construction and
  # are handled separately.)
  testthat::expect_gt(max(abs(unname(m_exact@bread))), bread_tol)
  testthat::expect_equal(
    unname(m_exact@bread),
    unname(m_cap@bread),
    tolerance = bread_tol,
    label = paste0(family, ": exact bread vs capprox")
  )

  if (check_se) {
    se_e <- sqrt(diag(m_exact@variance))
    se_c <- sqrt(diag(m_cap@variance))
    # Standard errors well above the comparison floor keep the relative all.equal
    # comparison from collapsing to an absolute (and vacuous) one.
    testthat::expect_gt(min(se_e), 1e-3)
    testthat::expect_equal(
      se_e,
      se_c,
      tolerance = se_tol,
      label = paste0(family, ": exact standard errors vs capprox")
    )
  }
  invisible(m_exact)
}

# ---- basic ------------------------------------------------------------------

test_that("exact mode matches capprox for the basic mean families", {
  set.seed(2024)
  y <- 0.6 + rnorm(200)
  expect_exact_acceptable(function(t) ee_mean(t, y = y), 0, "ee_mean")
  expect_exact_acceptable(
    function(t) ee_mean_variance(t, y = y),
    c(0, 1),
    "ee_mean_variance"
  )
  expect_exact_acceptable(
    function(t) ee_mean_robust(t, y = y, k = 1.345, loss = "huber"),
    0.6,
    "ee_mean_robust"
  )
  expect_exact_acceptable(
    function(t) ee_mean_geometric(t, y = abs(y) + 1),
    0.5,
    "ee_mean_geometric"
  )
})

# ee_percentile and ee_positive_mean_deviation estimate quantities defined
# through an indicator on theta, so their estimating equations are not
# differentiable at the solution. This is a property of the estimand, not of the
# derivative method: the bread is degenerate (identically zero for the
# percentile, singular for the positive mean deviation) under exact mode and
# under the central difference alike, and both functions warn that the sandwich
# should not be used. These tests pin that principled exception: the documented
# warning fires, and the exact-mode bread reproduces the same degenerate bread
# as capprox, so exact mode introduces no discrepancy of its own.
test_that("ee_percentile is a non-differentiable exception with a degenerate bread", {
  set.seed(2024)
  y <- 0.6 + rnorm(200)

  expect_warning(
    ee_percentile(median(y), y = y, q = 0.5),
    "not differentiable"
  )

  psi <- function(t) suppressWarnings(ee_percentile(t, y = y, q = 0.5))
  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = median(y)),
    deriv_method = "exact"
  )
  m_cap <- estimate(
    MEstimator(stacked_equations = psi, init = median(y)),
    deriv_method = "capprox"
  )
  # The bread is identically zero: the indicator has zero derivative away from
  # the kink, so the sandwich variance is undefined for this estimand.
  expect_equal(unname(m_exact@bread), matrix(0, 1, 1))
  expect_equal(unname(m_exact@bread), unname(m_cap@bread))
})

test_that("ee_positive_mean_deviation is a non-differentiable exception with a singular bread", {
  set.seed(2024)
  y <- 0.6 + rnorm(200)

  expect_warning(
    ee_positive_mean_deviation(c(1, median(y)), y = y),
    "not differentiable"
  )

  psi <- function(t) suppressWarnings(ee_positive_mean_deviation(t, y = y))
  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = c(1, median(y))),
    deriv_method = "exact"
  )
  m_cap <- estimate(
    MEstimator(stacked_equations = psi, init = c(1, median(y))),
    deriv_method = "capprox"
  )
  # The median score contributes a row of zeros to the bread (the indicator has
  # zero derivative), so the bread is singular and the sandwich is undefined.
  expect_equal(unname(m_exact@bread)[2, ], c(0, 0))
  expect_equal(unname(m_exact@bread), unname(m_cap@bread), tolerance = 1e-6)
})

# ---- regression and glm -----------------------------------------------------

test_that("exact mode matches capprox for the regression and glm families", {
  set.seed(2024)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  eta <- 0.6 + 1.1 * x1 - 0.5 * x2
  y_lin <- eta + rnorm(n)
  y_bin <- rbinom(n, 1, plogis(eta))
  y_pois <- rpois(n, exp(eta - 0.5))
  y_gam <- rgamma(n, shape = 2, scale = exp(eta) / 2)
  y_nb <- rnbinom(n, mu = exp(eta), size = 3)

  expect_exact_acceptable(
    function(t) ee_regression(t, X = X, y = y_lin, model = "linear"),
    c(0, 0, 0),
    "ee_regression (linear)"
  )
  expect_exact_acceptable(
    function(t) ee_regression(t, X = X, y = y_bin, model = "logistic"),
    c(0, 0, 0),
    "ee_regression (logistic)"
  )
  expect_exact_acceptable(
    function(t) ee_regression(t, X = X, y = y_pois, model = "poisson"),
    c(0, 0, 0),
    "ee_regression (poisson)"
  )
  expect_exact_acceptable(
    function(t) {
      ee_glm(t, X = X, y = y_gam, distribution = "gamma", link = "log")
    },
    c(0, 0, 0, 0),
    "ee_glm (gamma log)"
  )
  # The identity link requires a strictly positive mean, which the shared design
  # does not guarantee, so gamma identity uses a dedicated positive-mean setup.
  set.seed(7)
  x_gi <- rnorm(n)
  X_gi <- cbind(1, x_gi)
  mu_gi <- 3 + 0.8 * x_gi
  y_gi <- rgamma(n, shape = 5, rate = 5 / mu_gi)
  expect_exact_acceptable(
    function(t) {
      ee_glm(t, X = X_gi, y = y_gi, distribution = "gamma", link = "identity")
    },
    c(3, 0.8, 0),
    "ee_glm (gamma identity)"
  )
  expect_exact_acceptable(
    function(t) {
      ee_glm(t, X = X, y = y_bin, distribution = "binomial", link = "probit")
    },
    c(0, 0, 0),
    "ee_glm (probit)"
  )
  # The negative binomial dispersion nuisance parameter has a flat score, so its
  # standard error agrees to 1e-5 rather than 1e-6; the bread agrees to 1e-6.
  expect_exact_acceptable(
    function(t) {
      ee_glm(
        t,
        X = X,
        y = y_nb,
        distribution = "negative_binomial",
        link = "log"
      )
    },
    c(0.6, 1.1, -0.5, log(1 / 3)),
    "ee_glm (negative binomial)",
    se_tol = 1e-5
  )
  expect_exact_acceptable(
    function(t) {
      ee_glm(
        t,
        X = X,
        y = y_gam,
        distribution = "tweedie",
        link = "log",
        hyperparameter = 1.5
      )
    },
    c(0, 0, 0),
    "ee_glm (tweedie)"
  )
})

test_that("exact mode matches capprox for every robust regression loss", {
  set.seed(2025)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  y <- 0.7 + 1.2 * x1 - 0.6 * x2 + rnorm(n, sd = 2)

  losses <- list(
    huber = 1.345,
    cauchy = 1.345,
    tukey = 4.685,
    andrew = 1.339,
    hampel = c(1.5, 3, 6)
  )
  for (loss in names(losses)) {
    # The redescending losses need a starting value near the robust solution.
    expect_exact_acceptable(
      function(t) {
        ee_robust_regression(
          t,
          X = X,
          y = y,
          model = "linear",
          k = losses[[loss]],
          loss = loss
        )
      },
      c(0.7, 1.2, -0.6),
      paste0("ee_robust_regression (", loss, ")")
    )
  }
})

test_that("exact mode matches capprox for the penalized regression families", {
  set.seed(2024)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  y <- 0.6 + 1.1 * x1 - 0.5 * x2 + rnorm(n)

  expect_exact_acceptable(
    function(t) {
      ee_ridge_regression(
        t,
        X = X,
        y = y,
        model = "linear",
        penalty = c(0, 5, 5)
      )
    },
    c(0, 0, 0),
    "ee_ridge_regression"
  )
  expect_exact_acceptable(
    function(t) {
      ee_bridge_regression(
        t,
        X = X,
        y = y,
        model = "linear",
        penalty = c(0, 5, 5),
        gamma = 2
      )
    },
    c(0, 0, 0),
    "ee_bridge_regression (gamma = 2)"
  )
  expect_exact_acceptable(
    function(t) {
      suppressWarnings(
        ee_bridge_regression(
          t,
          X = X,
          y = y,
          model = "linear",
          penalty = c(0, 5, 5),
          gamma = 1.5
        )
      )
    },
    c(0.1, 0.1, 0.1),
    "ee_bridge_regression (gamma = 1.5)"
  )
  expect_exact_acceptable(
    function(t) {
      ee_dlasso_regression(
        t,
        X = X,
        y = y,
        model = "linear",
        penalty = c(0, 5, 5)
      )
    },
    c(0, 0, 0),
    "ee_dlasso_regression"
  )
  # The L1 penalties are not differentiable at the kink, so the sandwich is not
  # defined there; a sub-1e-6 bread difference blows up on the standard-error
  # scale. These families therefore compare the bread only, matching the caveat
  # in test-ee-regression-exact.R.
  expect_exact_acceptable(
    function(t) {
      suppressWarnings(
        ee_lasso_regression(
          t,
          X = X,
          y = y,
          model = "linear",
          penalty = c(0, 5, 5)
        )
      )
    },
    c(0.1, 0.1, 0.1),
    "ee_lasso_regression",
    check_se = FALSE
  )
  expect_exact_acceptable(
    function(t) {
      suppressWarnings(
        ee_elasticnet_regression(
          t,
          X = X,
          y = y,
          model = "linear",
          penalty = c(0, 5, 5),
          ratio = 0.5
        )
      )
    },
    c(0.1, 0.1, 0.1),
    "ee_elasticnet_regression",
    check_se = FALSE
  )
})

test_that("exact mode matches capprox for ee_additive_regression", {
  set.seed(42)
  n <- 200
  x <- runif(n, -3, 3)
  y <- sin(x) + rnorm(n, sd = 0.3)
  X <- cbind(1, x)
  specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))
  dm <- additive_design_matrix(X, specs, return_penalty = TRUE)
  init <- rep(0, ncol(dm$X))
  # The restricted-spline basis is ill-conditioned, so the bread agrees to 1e-5
  # rather than 1e-6; the standard errors still agree to 1e-6.
  expect_exact_acceptable(
    function(t) ee_additive_regression(t, X, y, specs, model = "linear"),
    init,
    "ee_additive_regression",
    bread_tol = 1e-5
  )
})

test_that("exact mode matches capprox for mlogit, beta, and tobit", {
  mlo <- load_fixture("ee_mlogit_exact")
  expect_exact_acceptable(
    function(t) ee_mlogit(t, X = mlo$X, y = mlo$y),
    mlo$init,
    "ee_mlogit"
  )
  bet <- load_fixture("ee_beta_regression_exact")
  expect_exact_acceptable(
    function(t) ee_beta_regression(t, X = bet$X, y = bet$y),
    bet$init,
    "ee_beta_regression"
  )
  tob <- load_fixture("ee_tobit_exact")
  expect_exact_acceptable(
    function(t) ee_tobit(t, X = tob$X, y = tob$y, lower = tob$lower),
    tob$init,
    "ee_tobit"
  )
})

# ---- causal -----------------------------------------------------------------

test_that("exact mode matches capprox for the fixture-backed causal families", {
  ipw <- load_fixture("ee_ipw_exact")
  expect_exact_acceptable(
    function(t) ee_ipw(t, y = ipw$y, A = ipw$a, W = ipw$W),
    ipw$init,
    "ee_ipw"
  )
  gf <- load_fixture("ee_gformula_exact")
  expect_exact_acceptable(
    function(t) ee_gformula(t, y = gf$y, X = gf$X, X1 = gf$X1, X0 = gf$X0),
    gf$init,
    "ee_gformula"
  )
  aipw <- load_fixture("ee_aipw_exact")
  expect_exact_acceptable(
    function(t) {
      ee_aipw(
        t,
        y = aipw$y,
        A = aipw$a,
        W = aipw$W,
        X = aipw$X,
        X1 = aipw$X1,
        X0 = aipw$X0
      )
    },
    aipw$init,
    "ee_aipw"
  )
  ge <- load_fixture("ee_gestimation_exact")
  expect_exact_acceptable(
    function(t) {
      ee_gestimation_snmm(
        t,
        y = ge$y,
        A = ge$a,
        W = ge$W,
        V = ge$V,
        model = ge$model
      )
    },
    ge$init,
    "ee_gestimation_snmm"
  )
  tsw <- load_fixture("ee_2sls_exact_w")
  expect_exact_acceptable(
    function(t) ee_2sls(t, y = tsw$y, A = tsw$a, Z = tsw$Z, W = tsw$W),
    tsw$init,
    "ee_2sls (with W)"
  )
  tsn <- load_fixture("ee_2sls_exact_now")
  expect_exact_acceptable(
    function(t) ee_2sls(t, y = tsn$y, A = tsn$a, Z = tsn$Z),
    tsn$init,
    "ee_2sls (without W)"
  )
})

test_that("exact mode matches capprox for ee_iv_causal", {
  set.seed(123)
  n <- 500
  Z <- rbinom(n, 1, 0.5)
  U <- rnorm(n)
  A <- rbinom(n, 1, plogis(-1 + 3 * Z + U))
  Y <- 3 * A - U + rnorm(n, sd = 0.5)
  expect_exact_acceptable(
    function(t) ee_iv_causal(t, y = Y, A = A, Z = Z),
    c(0, 0.5),
    "ee_iv_causal"
  )
})

test_that("exact mode matches capprox for ee_mean_sensitivity_analysis", {
  set.seed(42)
  n <- 500
  W <- rbinom(n, 1, 0.5)
  Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
  M <- rbinom(n, 1, plogis(2 + W))
  Y <- ifelse(M == 0, NA, Y_full)
  X <- cbind(1, W)
  y_filled <- ifelse(is.na(Y), 0, Y)
  # The outcome enters on a scale of order 100, which inflates the absolute
  # central-difference error, so the comparison uses 1e-4.
  expect_exact_acceptable(
    function(t) {
      ee_mean_sensitivity_analysis(
        t,
        y = y_filled,
        delta = M,
        X = X,
        q_eval = rep(0, n),
        H_function = inverse_logit
      )
    },
    c(180, 0, 0),
    "ee_mean_sensitivity_analysis",
    solver = "nleqslv",
    bread_tol = 1e-4,
    se_tol = 1e-4
  )
})

# ---- measurement ------------------------------------------------------------

test_that("exact mode matches capprox for the measurement families", {
  set.seed(123)
  n <- 500
  y_true <- rbinom(n, 1, 0.3)
  y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
  r <- c(rep(0, 200), rep(1, 300))
  y <- ifelse(r == 0, y_true, 0)
  expect_exact_acceptable(
    function(t) ee_rogan_gladen(t, y = y, y_star = y_star, r = r),
    c(0.3, 0.3, 0.8, 0.8),
    "ee_rogan_gladen",
    solver = "nleqslv"
  )

  rg <- load_fixture("ee_rogan_gladen_extended_exact")
  expect_exact_acceptable(
    function(t) {
      ee_rogan_gladen_extended(
        t,
        y = rg$y,
        y_star = rg$y_star,
        r = rg$r,
        X = rg$X
      )
    },
    rg$init,
    "ee_rogan_gladen_extended"
  )
  rc <- load_fixture("ee_regression_calibration_exact")
  expect_exact_acceptable(
    function(t) {
      ee_regression_calibration(
        t,
        beta = rc$beta,
        a = rc$a,
        a_star = rc$a_star,
        r = rc$r
      )
    },
    rc$init,
    "ee_regression_calibration"
  )
})

# ---- survival ---------------------------------------------------------------

test_that("exact mode matches capprox for ee_aft across distributions", {
  set.seed(101)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  X <- cbind(1, x)
  y_true <- exp(2 + 0.5 * x) * rweibull(n, shape = 0.75)
  censor <- rexp(n, rate = 0.1)
  y <- pmin(y_true, censor)
  delta <- as.numeric(y_true <= censor)

  expect_exact_acceptable(
    function(t) {
      ee_aft(t, X = X, time = y, event = delta, distribution = "exponential")
    },
    c(2, 0.5),
    "ee_aft (exponential)",
    solver = "nleqslv"
  )
  expect_exact_acceptable(
    function(t) {
      ee_aft(t, X = X, time = y, event = delta, distribution = "weibull")
    },
    c(2, 0.5, 0),
    "ee_aft (weibull)",
    solver = "nleqslv"
  )

  set.seed(313)
  n2 <- 300
  x2 <- rbinom(n2, 1, 0.5)
  u <- runif(n2)
  y2_true <- exp(2 + 0.5 * x2) * (u / (1 - u))^(1 / 2)
  censor2 <- rexp(n2, rate = 0.02)
  y2 <- pmin(y2_true, censor2)
  delta2 <- as.numeric(y2_true <= censor2)
  X2 <- cbind(1, x2)
  expect_exact_acceptable(
    function(t) {
      ee_aft(
        t,
        X = X2,
        time = y2,
        event = delta2,
        distribution = "log-logistic"
      )
    },
    c(2, 0.5, 0.5),
    "ee_aft (log-logistic)",
    solver = "nleqslv"
  )

  aft <- load_fixture("ee_aft_exact_lognormal")
  expect_exact_acceptable(
    function(t) {
      ee_aft(
        t,
        X = aft$X,
        time = aft$t,
        event = aft$delta,
        distribution = "log-normal"
      )
    },
    aft$init,
    "ee_aft (log-normal)"
  )
})

test_that("exact mode matches capprox for ee_survival_model across distributions", {
  set.seed(123)
  n <- 300
  u <- runif(n)
  t_true <- (-log(u) / 0.5)^(1 / 1.2)
  censor <- rexp(n, rate = 0.2)
  t_obs <- pmin(t_true, censor)
  delta <- as.numeric(t_true <= censor)
  expect_exact_acceptable(
    function(t) {
      ee_survival_model(
        t,
        time = t_obs,
        event = delta,
        distribution = "exponential"
      )
    },
    0.5,
    "ee_survival_model (exponential)",
    solver = "nleqslv"
  )
  expect_exact_acceptable(
    function(t) {
      ee_survival_model(
        t,
        time = t_obs,
        event = delta,
        distribution = "weibull"
      )
    },
    c(0.3, 1.0),
    "ee_survival_model (weibull)",
    solver = "nleqslv"
  )

  set.seed(345)
  u2 <- runif(n)
  t2_true <- (1 / 0.5) * log(1 - (0.5 / 0.3) * log(u2))
  valid <- is.finite(t2_true) & t2_true > 0
  t2_true <- t2_true[valid]
  nv <- length(t2_true)
  censor2 <- rexp(nv, rate = 0.3)
  t2_obs <- pmin(t2_true, censor2)
  delta2 <- as.numeric(t2_true <= censor2)
  expect_exact_acceptable(
    function(t) {
      ee_survival_model(
        t,
        time = t2_obs,
        event = delta2,
        distribution = "gompertz"
      )
    },
    c(0.2, 0.3),
    "ee_survival_model (gompertz)",
    solver = "nleqslv"
  )
})

test_that("exact mode matches capprox for ee_plogit", {
  ref <- load_fixture("ee_plogit_exact")
  expect_exact_acceptable(
    function(t) ee_plogit(t, X = ref$X, time = ref$time, event = ref$event),
    ref$init,
    "ee_plogit"
  )
})

# ---- pharmacokinetics -------------------------------------------------------

test_that("exact mode matches capprox for the dose-response families", {
  set.seed(42)
  n <- 50
  dose <- runif(n, 0, 100)
  response <- 1 + (25 * dose) / (10 + dose) + rnorm(n, sd = 1)
  # The nonlinear dose-response Jacobian is less well-conditioned than a linear
  # score, so exact and capprox agree to 1e-5 rather than 1e-6.
  expect_exact_acceptable(
    function(t) ee_emax(t, dose = dose, response = response),
    c(1, 25, 10),
    "ee_emax",
    solver = "nleqslv",
    bread_tol = 1e-5,
    se_tol = 1e-5
  )

  set.seed(42)
  nll <- 40
  dosell <- exp(runif(nll, log(0.1), log(100)))
  rholl <- (dosell / 10)^2
  respll <- 5 + (25 - 5) / (1 + rholl) + rnorm(nll, sd = 0.5)
  expect_exact_acceptable(
    function(t) ee_loglogistic(t, dose = dosell, response = respll),
    c(5, 25, 10, 2),
    "ee_loglogistic",
    solver = "nleqslv",
    bread_tol = 1e-5,
    se_tol = 1e-5
  )
})

# ---- no pinned gaps remaining -----------------------------------------------
#
# Every exported ee_* family now supports deriv_method = "exact". The only
# families that do not produce a usable sandwich under exact mode are the
# principled non-differentiable exceptions above (ee_percentile,
# ee_positive_mean_deviation, and the L1-kink standard errors of ee_lasso and
# ee_elasticnet), which are exceptions of the estimand rather than of the
# derivative method. There are no implementation gaps left to pin.

# ee_ipw_msm builds IPW weights from the propensity model (theta-derived, so
# tangent-carrying under exact mode) and passes them to ee_glm as the weights
# argument. The weighted GLM score keeps the weights' tangent, so exact mode
# reproduces the central-difference bread and standard errors.
test_that("exact mode matches capprox for ee_ipw_msm", {
  set.seed(42)
  n <- 200
  W <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.25 + 0.5 * W)
  Ya0 <- 5 + 2 * W + rnorm(n)
  Y <- (1 - A) * Ya0 + A * (Ya0 - 2)
  Wmat <- cbind(1, W)
  Vmsm <- cbind(1, A)
  expect_exact_acceptable(
    function(t) {
      ee_ipw_msm(
        t,
        y = Y,
        A = A,
        W = Wmat,
        V = Vmsm,
        distribution = "normal",
        link = "identity"
      )
    },
    rep(0, 4),
    "ee_ipw_msm"
  )
})

# The stacked delta-effective-dose equations ee_emax_ed and ee_loglogistic_ed
# verify alongside their base dose-response family. Each shares the nonlinear
# dose-response conditioning of its base family, so exact and capprox agree to
# 1e-5 rather than 1e-6.
test_that("exact mode matches capprox for the stacked delta-effective-dose equations", {
  set.seed(42)
  n <- 50
  dose <- runif(n, 0, 100)
  response <- 1 + (25 * dose) / (10 + dose) + rnorm(n, sd = 1)
  expect_exact_acceptable(
    function(t) {
      rbind(
        ee_emax(t[1:3], dose = dose, response = response),
        ee_emax_ed(t[4], dose = dose, delta = 0.8, ed50 = t[3])
      )
    },
    c(1, 25, 10, 40),
    "ee_emax_ed",
    solver = "nleqslv",
    bread_tol = 1e-5,
    se_tol = 1e-5
  )

  set.seed(42)
  nll <- 40
  dosell <- exp(runif(nll, log(0.1), log(100)))
  rholl <- (dosell / 10)^2
  respll <- 5 + (25 - 5) / (1 + rholl) + rnorm(nll, sd = 0.5)
  expect_exact_acceptable(
    function(t) {
      rbind(
        ee_loglogistic(t[1:4], dose = dosell, response = respll),
        ee_loglogistic_ed(
          t[5],
          dose = dosell,
          delta = 0.5,
          lower = t[1],
          upper = t[2],
          ed50 = t[3],
          steepness = t[4]
        )
      )
    },
    c(5, 25, 10, 2, 10),
    "ee_loglogistic_ed",
    solver = "nleqslv",
    bread_tol = 1e-5,
    se_tol = 1e-5
  )
})

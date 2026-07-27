# Comprehensive cross-validation against Python Delicatessen fixtures
#
# Each test validates theta, full variance matrix, asymptotic variance,
# bread matrix, and confidence intervals against pre-computed Python results.

# ---- basic estimating equations ----

test_that("ee_mean matches Python Delicatessen", {
  ref <- load_fixture("ee_mean")
  y <- ref$y

  psi <- function(theta) {
    ee_mean(theta, y = y)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mean", tolerance = 1e-6)
})

test_that("ee_mean_variance matches Python Delicatessen", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    ee_mean_variance(theta, y = y)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mean_variance", tolerance = 1e-6)
})

test_that("ee_mean_robust (Huber) matches Python Delicatessen", {
  ref <- load_fixture("ee_mean_robust_huber")
  y <- ref$y
  k <- ref$k

  psi <- function(theta) {
    ee_mean_robust(theta, y = y, k = k, loss = "huber")
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mean_robust_huber", tolerance = 1e-6)
})

test_that("ee_mean_geometric matches Python Delicatessen", {
  ref <- load_fixture("ee_mean_geometric")
  y <- ref$y

  psi <- function(theta) {
    ee_mean_geometric(theta, y = y)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mean_geometric", tolerance = 1e-6)
})

test_that("ee_positive_mean_deviation matches Python Delicatessen", {
  ref <- load_fixture("ee_positive_mean_deviation")
  y <- ref$y

  psi <- function(theta) {
    suppressWarnings(ee_positive_mean_deviation(theta, y = y))
  }

  # The median equation uses an indicator function (non-differentiable),
  # so the EE has no exact root and solvers find different approximate
  # solutions depending on init and algorithm. Instead of comparing
  # solved theta, verify the EE formulation and sandwich components
  # match Python when evaluated at the same point.
  m <- MEstimator(stacked_equations = psi, init = ref$theta)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)
  expect_equal(unname(m@variance), as.matrix(ref$variance), tolerance = 1e-4)
  expect_equal(unname(m@bread), as.matrix(ref$bread), tolerance = 1e-4)
})

test_that("ee_percentile matches Python Delicatessen", {
  ref <- load_fixture("ee_percentile_50")
  y <- ref$y
  q <- ref$q

  psi <- function(theta) {
    suppressWarnings(ee_percentile(theta, y = y, q = q))
  }

  # Percentile EE is non-differentiable (indicator function).
  # Python uses solver='hybr', tolerance=1e-3, dx=1, deriv_method='fapprox'
  # to get approximate values. We match those settings here.
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(
    m,
    solver = "nleqslv",
    tolerance = 1e-3,
    dx = 1,
    deriv_method = "fapprox"
  )

  expect_python_match(m, "ee_percentile_50", tolerance = 1e-3)
})

# ---- causal ----

test_that("ee_gformula matches Python Delicatessen", {
  ref <- load_fixture("ee_gformula")
  X <- ref$X
  X1 <- ref$X1
  X0 <- ref$X0
  y <- ref$y

  psi <- function(theta) {
    ee_gformula(theta, y = y, X = X, X1 = X1, X0 = X0)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_gformula", tolerance = 1e-5)
})

test_that("ee_ipw matches Python Delicatessen", {
  ref <- load_fixture("ee_ipw")
  y <- ref$y
  A <- ref$a
  W <- ref$W

  psi <- function(theta) {
    ee_ipw(theta, y = y, A = A, W = W)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_ipw", tolerance = 1e-5)
})

test_that("ee_aipw matches Python Delicatessen", {
  ref <- load_fixture("ee_aipw")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  X <- ref$X
  X1 <- ref$X1
  X0 <- ref$X0

  psi <- function(theta) {
    ee_aipw(theta, y = y, A = A, W = W, X = X, X1 = X1, X0 = X0)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_aipw", tolerance = 1e-5)
})

test_that("ee_ipw_msm matches Python Delicatessen", {
  ref <- load_fixture("ee_ipw_msm")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  V <- ref$V

  psi <- function(theta) {
    ee_ipw_msm(
      theta,
      y = y,
      A = A,
      W = W,
      V = V,
      distribution = ref$distribution,
      link = ref$link
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_ipw_msm", tolerance = 1e-5)
})

test_that("ee_gestimation_snmm matches Python Delicatessen", {
  ref <- load_fixture("ee_gestimation_snmm")
  y <- ref$y
  A <- ref$a
  W <- ref$W
  V <- ref$V

  psi <- function(theta) {
    ee_gestimation_snmm(theta, y = y, A = A, W = W, V = V, model = ref$model)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_gestimation_snmm", tolerance = 1e-5)
})

test_that("ee_iv_causal matches Python Delicatessen", {
  ref <- load_fixture("ee_iv_causal")
  y <- ref$y
  A <- ref$a
  Z <- ref$z

  psi <- function(theta) {
    ee_iv_causal(theta, y = y, A = A, Z = Z)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_iv_causal", tolerance = 1e-5)
})

test_that("ee_2sls matches Python Delicatessen", {
  ref <- load_fixture("ee_2sls")
  y <- ref$y
  A <- ref$a
  Z <- ref$z

  psi <- function(theta) {
    ee_2sls(theta, y = y, A = A, Z = Z)
  }

  # 2SLS needs non-zero init to avoid singular first-stage predictions.
  # Python's LM solver handles [0,0] but R solvers don't.
  init <- c(1, 0.5)
  m <- MEstimator(stacked_equations = psi, init = init)
  m <- estimate(m)

  expect_python_match(m, "ee_2sls", tolerance = 1e-5)
})

test_that("ee_mean_sensitivity_analysis matches Python Delicatessen", {
  ref <- load_fixture("ee_mean_sensitivity_analysis")
  y <- ref$y
  delta <- ref$delta
  X <- ref$X
  q_eval <- ref$q_eval

  psi <- function(theta) {
    ee_mean_sensitivity_analysis(
      theta,
      y = y,
      delta = delta,
      X = X,
      q_eval = q_eval,
      H_function = inverse_logit
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mean_sensitivity_analysis", tolerance = 1e-4)
})

# ---- regression ----

test_that("ee_regression linear matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_linear")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_regression_linear", tolerance = 1e-5)
})

test_that("ee_regression logistic matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_logistic")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "logistic")
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_regression_logistic", tolerance = 1e-5)
})

test_that("ee_regression poisson matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_poisson")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "poisson")
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_regression_poisson", tolerance = 1e-5)
})

test_that("ee_ridge_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_ridge_regression")
  X <- ref$X
  y <- ref$y
  penalty <- ref$penalty

  psi <- function(theta) {
    ee_ridge_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      penalty = penalty
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_ridge_regression", tolerance = 1e-5)
})

test_that("ee_robust_regression (Huber) matches Python Delicatessen", {
  ref <- load_fixture("ee_robust_regression_huber")
  X <- ref$X
  y <- ref$y
  k <- ref$k

  psi <- function(theta) {
    ee_robust_regression(
      theta,
      X = X,
      y = y,
      model = "linear",
      k = k,
      loss = "huber"
    )
  }

  # Huber regression clips residuals beyond k, so from a starting point where
  # every residual is clipped the Jacobian is singular and Newton solvers
  # stall. Python delicatessen solves with scipy's Levenberg-Marquardt; the
  # "lm" solver wraps the same MINPACK routine and converges from the origin.
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "lm")

  expect_python_match(m, "ee_robust_regression_huber", tolerance = 1e-5)
})

# ---- penalized regression ----

test_that("ee_lasso_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_lasso_regression")
  X <- ref$X
  y <- ref$y

  # The L1 penalty warns on every evaluation (see test-ee-penalized-regression.R)
  psi <- function(theta) {
    suppressWarnings(
      ee_lasso_regression(
        theta,
        X = X,
        y = y,
        model = ref$model,
        penalty = ref$penalty
      )
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- suppressWarnings(estimate(m))

  expect_python_match(m, "ee_lasso_regression", tolerance = 1e-5)
})

test_that("ee_dlasso_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_dlasso_regression")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_dlasso_regression(
      theta,
      X = X,
      y = y,
      model = ref$model,
      penalty = ref$penalty
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_dlasso_regression", tolerance = 1e-5)
})

test_that("ee_elasticnet_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_elasticnet_regression")
  X <- ref$X
  y <- ref$y

  # The L1 component warns on every evaluation
  psi <- function(theta) {
    suppressWarnings(
      ee_elasticnet_regression(
        theta,
        X = X,
        y = y,
        model = ref$model,
        penalty = ref$penalty,
        ratio = ref$ratio
      )
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- suppressWarnings(estimate(m))

  expect_python_match(m, "ee_elasticnet_regression", tolerance = 1e-5)
})

test_that("ee_bridge_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_bridge_regression")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_bridge_regression(
      theta,
      X = X,
      y = y,
      model = ref$model,
      penalty = ref$penalty,
      gamma = ref$gamma
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_bridge_regression", tolerance = 1e-5)
})

# ---- GLM ----

test_that("ee_glm normal/identity matches Python Delicatessen", {
  ref <- load_fixture("ee_glm_normal_identity")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_glm(
      theta,
      X = X,
      y = y,
      distribution = ref$distribution,
      link = ref$link
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_glm_normal_identity", tolerance = 1e-5)
})

test_that("ee_glm poisson/log matches Python Delicatessen", {
  ref <- load_fixture("ee_glm_poisson_log")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_glm(
      theta,
      X = X,
      y = y,
      distribution = ref$distribution,
      link = ref$link
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_glm_poisson_log", tolerance = 1e-5)
})

test_that("ee_mlogit matches Python Delicatessen", {
  ref <- load_fixture("ee_mlogit")
  X <- ref$X
  y <- ref$y # one-hot encoded matrix

  psi <- function(theta) {
    ee_mlogit(theta, X = X, y = y)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_mlogit", tolerance = 1e-5)
})

test_that("ee_beta_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_beta_regression")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_beta_regression(theta, X = X, y = y)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_beta_regression", tolerance = 1e-5)
})

test_that("ee_tobit matches Python Delicatessen", {
  ref <- load_fixture("ee_tobit")
  X <- ref$X
  y <- ref$y

  psi <- function(theta) {
    ee_tobit(theta, X = X, y = y, lower = ref$lower)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_tobit", tolerance = 1e-5)
})

test_that("ee_additive_regression matches Python Delicatessen", {
  ref <- load_fixture("ee_additive_regression")
  X <- ref$X
  y <- ref$y

  # Reconstruct specifications from fixture because jsonlite can mangle
  # nested list structures when simplifyVector = TRUE
  specs <- list(
    NULL,
    list(
      knots = c(-2, -1, 0, 1, 2),
      penalty = 3,
      natural = TRUE,
      power = 3,
      normalized = FALSE
    )
  )

  psi <- function(theta) {
    ee_additive_regression(
      theta,
      X = X,
      y = y,
      specifications = specs,
      model = ref$model
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_additive_regression", tolerance = 1e-4)
})

# ---- survival ----

test_that("ee_survival_model Weibull matches Python Delicatessen", {
  ref <- load_fixture("ee_survival_model_weibull")
  t_obs <- ref$t
  delta <- ref$delta

  psi <- function(theta) {
    ee_survival_model(
      theta,
      time = t_obs,
      event = delta,
      distribution = ref$distribution
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_survival_model_weibull", tolerance = 1e-5)
})

test_that("ee_aft Weibull matches Python Delicatessen", {
  ref <- load_fixture("ee_aft_weibull")
  X <- ref$X
  t_obs <- ref$t
  delta <- ref$delta

  psi <- function(theta) {
    ee_aft(
      theta,
      X = X,
      time = t_obs,
      event = delta,
      distribution = ref$distribution
    )
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_aft_weibull", tolerance = 1e-4)
})

test_that("ee_plogit matches Python Delicatessen", {
  ref <- load_fixture("ee_plogit")
  X <- ref$X
  # In R ee_plogit: time = observed times, event = event indicator
  y_times <- ref$t
  t_delta <- ref$delta
  S <- ref$S

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y_times, event = t_delta, S = S)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_plogit", tolerance = 1e-4)
})

# ---- measurement ----

test_that("ee_rogan_gladen matches Python Delicatessen", {
  ref <- load_fixture("ee_rogan_gladen")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r

  psi <- function(theta) {
    ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_rogan_gladen", tolerance = 1e-5)
})

test_that("ee_rogan_gladen_extended matches Python Delicatessen", {
  ref <- load_fixture("ee_rogan_gladen_extended")
  y <- ref$y
  y_star <- ref$y_star
  r <- ref$r
  X <- ref$X

  psi <- function(theta) {
    ee_rogan_gladen_extended(theta, y = y, y_star = y_star, r = r, X = X)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_python_match(m, "ee_rogan_gladen_extended", tolerance = 1e-5)
})

test_that("ee_regression_calibration matches Python Delicatessen", {
  ref <- load_fixture("ee_regression_calibration")
  a <- ref$a
  a_star <- ref$a_star
  y <- ref$y
  r <- ref$r

  # Stacked: calibration (3 params) + naive logistic (2 params)
  X_main <- cbind(1, a_star)

  psi <- function(theta) {
    theta_calib <- theta[1:3]
    theta_main <- theta[4:5]

    # Naive logistic (main study only)
    y_safe <- ifelse(r == 0, -999, y)
    ee_logit <- ee_regression(
      theta_main,
      X = X_main,
      y = y_safe,
      model = "logistic"
    )
    ee_logit <- ee_logit * rep(r, each = nrow(ee_logit))

    # Regression calibration
    ee_calib <- ee_regression_calibration(
      theta_calib,
      beta = theta_main[2],
      a = a,
      a_star = a_star,
      r = r
    )

    rbind(ee_calib, ee_logit)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_regression_calibration", tolerance = 1e-5)
})

# ---- pharmacokinetics ----

test_that("ee_emax matches Python Delicatessen", {
  ref <- load_fixture("ee_emax")
  dose <- ref$dose
  response <- ref$response

  psi <- function(theta) {
    ee_emax(theta, dose = dose, response = response)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_emax", tolerance = 1e-4)
})

test_that("ee_loglogistic matches Python Delicatessen", {
  ref <- load_fixture("ee_loglogistic")
  dose <- ref$dose
  response <- ref$response

  psi <- function(theta) {
    ee_loglogistic(theta, dose = dose, response = response)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_loglogistic", tolerance = 1e-4)
})

test_that("ee_emax + ee_emax_ed stacked matches Python Delicatessen", {
  ref <- load_fixture("ee_emax_ed")
  dose <- ref$dose
  response <- ref$response
  delta <- ref$delta

  psi <- function(theta) {
    emax_ee <- ee_emax(theta[1:3], dose = dose, response = response)
    ed_ee <- ee_emax_ed(theta[4], dose = dose, delta = delta, ed50 = theta[3])
    rbind(emax_ee, ed_ee)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_emax_ed", tolerance = 1e-4)
})

test_that("ee_loglogistic + ee_loglogistic_ed stacked matches Python Delicatessen", {
  ref <- load_fixture("ee_loglogistic_ed")
  dose <- ref$dose
  response <- ref$response
  lower <- ref$lower
  delta <- ref$delta

  psi <- function(theta) {
    # 3PL: fix lower, estimate upper (theta[1]), ed50 (theta[2]),
    # steepness (theta[3])
    ll_ee <- ee_loglogistic(
      c(lower, theta[1:3]),
      dose = dose,
      response = response
    )
    ll_ee <- ll_ee[-1, , drop = FALSE] # Drop lower limit EE

    ed_ee <- ee_loglogistic_ed(
      theta[4],
      dose = dose,
      delta = delta,
      lower = lower,
      upper = theta[1],
      ed50 = theta[2],
      steepness = theta[3]
    )
    rbind(ll_ee, ed_ee)
  }

  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, solver = "nleqslv")

  expect_python_match(m, "ee_loglogistic_ed", tolerance = 1e-4)
})

# Tests for regression_predictions

test_that("regression_predictions returns correct shape", {
  set.seed(42)
  n_pred <- 10
  X <- cbind(1, rnorm(n_pred))
  theta <- c(1, 0.5)
  cov_mat <- matrix(c(0.1, 0.01, 0.01, 0.05), nrow = 2)

  result <- regression_predictions(X, theta, cov_mat)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), n_pred)
  expect_equal(ncol(result), 4)
  expect_named(result, c("predicted", "variance", "lower", "upper"))
})

test_that("regression_predictions computes correct predictions", {
  X <- cbind(1, c(0, 1, 2))
  theta <- c(2, 3)
  cov_mat <- diag(2) * 0.01

  result <- regression_predictions(X, theta, cov_mat)

  # y = 2 + 3*x -> 2, 5, 8
  expect_equal(result$predicted, c(2, 5, 8))

  # Variance: diag(X %*% cov %*% t(X))
  expected_var <- rowSums((X %*% cov_mat) * X)
  expect_equal(result$variance, expected_var)
})

test_that("regression_predictions CI contains predicted value", {
  set.seed(123)
  X <- cbind(1, rnorm(5))
  theta <- c(1, 0.5)
  cov_mat <- diag(2) * 0.1

  result <- regression_predictions(X, theta, cov_mat)
  expect_true(all(result$lower <= result$predicted))
  expect_true(all(result$upper >= result$predicted))
})

test_that("regression_predictions wider CI with larger alpha", {
  X <- cbind(1, c(0, 1))
  theta <- c(1, 1)
  cov_mat <- diag(2) * 0.1

  ci_95 <- regression_predictions(X, theta, cov_mat, alpha = 0.05)
  ci_99 <- regression_predictions(X, theta, cov_mat, alpha = 0.01)

  # 99% CI should be wider than 95% CI
  width_95 <- ci_95$upper - ci_95$lower
  width_99 <- ci_99$upper - ci_99$lower
  expect_true(all(width_99 > width_95))
})

test_that("regression_predictions supports offset", {
  X <- cbind(1, c(0, 1))
  theta <- c(1, 2)
  cov_mat <- diag(2) * 0.01
  off <- c(0.5, 0.5)

  result_no_off <- regression_predictions(X, theta, cov_mat)
  result_off <- regression_predictions(X, theta, cov_mat, offset = off)

  expect_equal(result_off$predicted, result_no_off$predicted + 0.5)
})

test_that("regression_predictions errors on invalid alpha", {
  X <- cbind(1, 1:3)
  theta <- c(0, 0)
  cov_mat <- diag(2)

  expect_error(regression_predictions(X, theta, cov_mat, alpha = 0), "alpha")
  expect_error(regression_predictions(X, theta, cov_mat, alpha = 1), "alpha")
})

test_that("regression_predictions works with MEstimator output", {
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  X <- cbind(1, x)

  psi <- function(theta) {
    ee_regression(theta, X = X, y = y, model = "linear")
  }
  m <- MEstimator(stacked_equations = psi, init = c(0, 0))
  m <- estimate(m)

  # Predict at new values
  x_new <- seq(-2, 2, length.out = 5)
  X_new <- cbind(1, x_new)
  preds <- regression_predictions(X_new, m@theta, m@variance)

  expect_equal(nrow(preds), 5)
  # Predictions should be close to true line
  expect_equal(preds$predicted, 1 + 2 * x_new, tolerance = 0.5)
  # All variances positive

  expect_true(all(preds$variance > 0))
})

# convert_survival_measures tests -------------------------------------------

test_that("convert_survival_measures returns correct survival", {
  expect_equal(convert_survival_measures(0.8, 0.01, "survival"), 0.8)
})

test_that("convert_survival_measures returns correct risk", {
  expect_equal(convert_survival_measures(0.8, 0.01, "risk"), 0.2)
})

test_that("convert_survival_measures returns correct cumulative hazard", {
  expect_equal(convert_survival_measures(0.8, 0.01, "cumulative_hazard"),
               -log(0.8))
})

test_that("convert_survival_measures returns correct density", {
  expect_equal(convert_survival_measures(0.8, 0.01, "density"), 0.01 * 0.8)
})

test_that("convert_survival_measures errors on invalid measure", {
  expect_error(convert_survival_measures(0.8, 0.01, "invalid"), "not supported")
})

test_that("convert_survival_measures works with vectors", {
  surv <- c(0.9, 0.8, 0.7)
  haz <- c(0.01, 0.02, 0.03)
  result <- convert_survival_measures(surv, haz, "risk")
  expect_equal(result, 1 - surv)
})

# survival_predictions tests -----------------------------------------------

test_that("survival_predictions returns correct shape", {
  # Use known Weibull parameters
  theta <- c(0.5, 1.2)  # lambda, gamma
  cov_mat <- diag(2) * 0.01
  times <- c(0.5, 1, 1.5, 2)

  result <- survival_predictions(times, theta, cov_mat,
                                  distribution = "weibull")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 4)
  expect_named(result, c("time", "predicted", "variance", "lower", "upper"))
})

test_that("survival_predictions survival is between 0 and 1", {
  theta <- c(0.5, 1.5)
  cov_mat <- diag(2) * 0.001
  times <- seq(0.1, 3, by = 0.5)

  result <- survival_predictions(times, theta, cov_mat,
                                  distribution = "weibull",
                                  measure = "survival")
  expect_true(all(result$predicted >= 0 & result$predicted <= 1))
})

test_that("survival_predictions with ee_survival_model", {
  set.seed(42)
  n <- 200
  # Generate Weibull data: S(t) = exp(-lambda * t^gamma)
  lambda_true <- 0.5
  gamma_true <- 1.2
  u <- runif(n)
  t_event <- (-log(u) / lambda_true)^(1 / gamma_true)
  # Simple right censoring
  t_censor <- rexp(n, rate = 0.2)
  t_obs <- pmin(t_event, t_censor)
  delta <- as.numeric(t_event <= t_censor)

  psi <- function(theta) {
    ee_survival_model(theta, time = t_obs, event = delta,
                       distribution = "weibull")
  }
  m <- MEstimator(stacked_equations = psi, init = c(0.5, 1))
  m <- estimate(m, solver = "nleqslv")

  # Predict survival at several time points
  pred_times <- c(0.5, 1, 2)
  preds <- survival_predictions(pred_times, m@theta, m@variance,
                                 distribution = "weibull")

  expect_equal(nrow(preds), 3)
  expect_true(all(preds$variance > 0))
  # Survival should decrease over time
  expect_true(preds$predicted[1] > preds$predicted[3])
})

test_that("survival_predictions errors on invalid alpha", {
  expect_error(
    survival_predictions(1, c(0.5, 1), diag(2), "weibull", alpha = 0),
    "alpha"
  )
})

# survival_predictions exact-derivative tests (bd-1v9o.9) -------------------
#
# Python's survival_predictions computes the delta-method variance with exact
# forward-mode automatic differentiation: it calls delta_method, which defaults
# to deriv_method='exact'. The R port takes a `deriv_method` argument that it
# forwards to delta_method. Its default is "capprox" (central difference), so
# the no-argument behavior is central-difference numerical differentiation;
# deriv_method = "exact" selects the autodiff path that Python uses internally.
#
# The delta-method g_func concatenates the per-time measures with c() rather
# than coercing with vapply(..., numeric(1)): under exact autodiff the per-time
# values are tangent-carrying pair objects, not plain doubles, so a
# vapply(numeric(1)) template cannot hold them. The point-estimate line still
# uses vapply, which is correct because it evaluates at the numeric theta where
# no tangents are present.
#
# Fixture provenance: fixtures/survival_predictions.json, produced by
# fixtures/generate_survival_predictions_fixture.py. That script fits
# parametric survival models with ee_survival_model on the deterministic
# Collett (2015) breast cancer data (no RNG) for the weibull and exponential
# distributions, records theta and covariance, then records
# survival_predictions output (predicted measure, delta-method variance, and
# Wald confidence limits) across five measures, five times, and two alpha
# levels. Because Python computes the variance with exact autodiff, the
# recorded variances and confidence limits pin the exact-mode values. Feeding
# the recorded theta and covariance to R isolates the prediction and
# delta-method logic from survival-model fitting parity.
#
# Tolerances: point predictions are analytic and do not depend on the
# derivative method; they are pinned at 1e-6 (they agree with Python to
# machine precision, ~0 relative on this fixture). For deriv_method = "exact"
# the variance is compared on the SE (square-root) scale and the confidence
# limits directly, both at 1e-8: R and Python both evaluate the same
# closed-form Jacobian with forward-mode autodiff, so they agree far tighter
# than 1e-8 (rounding only). For deriv_method = "capprox" the same quantities
# are pinned at 1e-5, matching the aft_predictions_function convention; R's
# central-difference SE agrees with Python's exact SE to ~2e-7 relative on
# this fixture, so 1e-5 leaves margin for that documented method difference.
# The SE scale keeps the variance comparison in all.equal's relative mode:
# several measures have variances as small as ~3e-9 (density) that would fall
# into absolute mode (where a variance of 0 passes vacuously) if compared
# directly, whereas every SE vector has mean magnitude of at least ~4e-4.

test_that("survival_predictions with deriv_method = 'exact' matches Python", {
  ref <- load_fixture("survival_predictions")

  for (i in seq_len(nrow(ref$scenarios))) {
    s <- ref$scenarios[i, ]
    dist <- s$distribution
    measure <- s$measure
    theta <- ref$models[[dist]]$theta
    covariance <- as.matrix(ref$models[[dist]]$covariance)
    times <- ref$scenarios$times[[i]]
    alpha <- s$alpha

    result <- survival_predictions(
      times = times, theta = theta, covariance = covariance,
      distribution = dist, measure = measure, alpha = alpha,
      deriv_method = "exact"
    )

    lbl <- paste0(dist, "/", measure, "/alpha=", alpha)

    expect_equal(
      result$time, times,
      tolerance = 1e-6, label = paste0(lbl, ": time")
    )
    expect_equal(
      result$predicted, ref$scenarios$predicted[[i]],
      tolerance = 1e-6, label = paste0(lbl, ": predicted")
    )
    # Variance compared on the SE scale to stay in relative mode; see the
    # section header for why direct comparison would be vacuous for the
    # small-variance measures.
    expect_equal(
      sqrt(result$variance), sqrt(ref$scenarios$variance[[i]]),
      tolerance = 1e-8, label = paste0(lbl, ": se")
    )
    expect_equal(
      result$lower, ref$scenarios$lower[[i]],
      tolerance = 1e-8, label = paste0(lbl, ": lower")
    )
    expect_equal(
      result$upper, ref$scenarios$upper[[i]],
      tolerance = 1e-8, label = paste0(lbl, ": upper")
    )
  }
})

test_that("survival_predictions exact-mode variances are positive without clamping", {
  # survival_predictions forms its SE with sqrt(pmax(variance, 0)), a safeguard
  # for numerical differentiation. Exact autodiff makes the safeguard
  # unnecessary: the delta-method variance is diag(G Sigma G^T) with Sigma a
  # valid (positive semi-definite) covariance, so it is nonnegative by
  # construction, and it is strictly positive for every scenario in the
  # fixture. Asserting strict positivity pins that the clamp is a no-op here
  # (pmax(v, 0) == v), so the reported values come from the Jacobian and not
  # from the clamp; matching the strictly positive Python fixture keeps this
  # from passing vacuously (a clamp-produced 0 would fail both checks).
  ref <- load_fixture("survival_predictions")

  for (i in seq_len(nrow(ref$scenarios))) {
    s <- ref$scenarios[i, ]
    dist <- s$distribution
    measure <- s$measure
    theta <- ref$models[[dist]]$theta
    covariance <- as.matrix(ref$models[[dist]]$covariance)

    result <- survival_predictions(
      times = ref$scenarios$times[[i]], theta = theta,
      covariance = covariance, distribution = dist,
      measure = measure, alpha = s$alpha, deriv_method = "exact"
    )

    lbl <- paste0(dist, "/", measure, "/alpha=", s$alpha)
    expect_true(
      all(result$variance > 0),
      label = paste0(lbl, ": all variances strictly positive")
    )
    # The Python fixture variances are all strictly positive; matching them on
    # the SE scale confirms the R exact-mode variance equals Python's, so the
    # positivity is genuine rather than a clamp artifact.
    expect_equal(
      sqrt(result$variance), sqrt(ref$scenarios$variance[[i]]),
      tolerance = 1e-8, label = paste0(lbl, ": se matches positive Python")
    )
  }
})

test_that("survival_predictions with deriv_method = 'capprox' matches Python", {
  # The capprox path is the survival_predictions default and must keep matching
  # Python (at the looser numerical-vs-exact tolerance) once the deriv_method
  # argument is wired through, so both derivative modes are pinned.
  ref <- load_fixture("survival_predictions")

  for (i in seq_len(nrow(ref$scenarios))) {
    s <- ref$scenarios[i, ]
    dist <- s$distribution
    measure <- s$measure
    theta <- ref$models[[dist]]$theta
    covariance <- as.matrix(ref$models[[dist]]$covariance)
    times <- ref$scenarios$times[[i]]

    result <- survival_predictions(
      times = times, theta = theta, covariance = covariance,
      distribution = dist, measure = measure, alpha = s$alpha,
      deriv_method = "capprox"
    )

    lbl <- paste0(dist, "/", measure, "/alpha=", s$alpha)
    expect_equal(
      result$predicted, ref$scenarios$predicted[[i]],
      tolerance = 1e-6, label = paste0(lbl, ": predicted")
    )
    expect_equal(
      sqrt(result$variance), sqrt(ref$scenarios$variance[[i]]),
      tolerance = 1e-5, label = paste0(lbl, ": se")
    )
    expect_equal(
      result$lower, ref$scenarios$lower[[i]],
      tolerance = 1e-5, label = paste0(lbl, ": lower")
    )
    expect_equal(
      result$upper, ref$scenarios$upper[[i]],
      tolerance = 1e-5, label = paste0(lbl, ": upper")
    )
  }
})

test_that("survival_predictions default derivative method reproduces capprox", {
  # Pins that the deriv_method default stays "capprox" so existing behavior is
  # unchanged: the no-argument call must equal an explicit deriv_method =
  # "capprox" call. (The explicit call is the target contract; this test starts
  # passing once bd-1v9o.10 adds the argument.)
  ref <- load_fixture("survival_predictions")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)
  times <- ref$times

  default_call <- survival_predictions(
    times = times, theta = theta, covariance = covariance,
    distribution = "weibull", measure = "risk"
  )
  capprox_call <- survival_predictions(
    times = times, theta = theta, covariance = covariance,
    distribution = "weibull", measure = "risk", deriv_method = "capprox"
  )

  expect_equal(default_call, capprox_call)
})

# aft_predictions_individual tests -----------------------------------------

test_that("aft_predictions_individual returns correct shape", {
  set.seed(42)
  n <- 20
  X <- cbind(1, rnorm(n))
  theta <- c(1, 0.5, log(1))  # beta0, beta1, log(1/sigma)
  times <- c(0.5, 1, 2)

  result <- aft_predictions_individual(X, times, theta,
                                        distribution = "weibull")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), n)
  expect_equal(ncol(result), 3)
})

test_that("aft_predictions_individual survival is in [0, 1]", {
  set.seed(42)
  n <- 10
  X <- cbind(1, rnorm(n))
  theta <- c(2, 0.3, log(1.5))
  times <- c(1, 5, 10)

  result <- aft_predictions_individual(X, times, theta,
                                        distribution = "weibull",
                                        measure = "survival")
  expect_true(all(unlist(result) >= 0))
  expect_true(all(unlist(result) <= 1))
})

test_that("aft_predictions_individual survival decreases over time", {
  X <- cbind(1, 0)  # single individual
  theta <- c(2, 0.3, log(1.5))
  times <- c(1, 5, 10, 20)

  result <- aft_predictions_individual(X, times, theta,
                                        distribution = "weibull",
                                        measure = "survival")
  surv <- as.numeric(result[1, ])
  # Survival should be non-increasing over time
  expect_true(all(diff(surv) <= 0))
})

test_that("aft_predictions_individual works with MEstimator", {
  set.seed(123)
  n <- 200
  X <- cbind(1, rbinom(n, 1, 0.5))
  beta_true <- c(2, 0.5)
  sigma_true <- 0.8
  eps <- rnorm(n)
  t_event <- exp(as.numeric(X %*% beta_true) + sigma_true * eps)
  t_censor <- rexp(n, rate = 0.05)
  t_obs <- pmin(t_event, t_censor)
  delta <- as.numeric(t_event <= t_censor)

  psi <- function(theta) {
    ee_aft(theta, X = X, time = t_obs, event = delta,
           distribution = "log-normal")
  }
  m <- MEstimator(stacked_equations = psi,
                   init = c(mean(log(t_obs)), 0, log(1 / sd(log(t_obs)))))
  m <- estimate(m, solver = "nleqslv")

  pred_times <- c(1, 5, 10)
  preds <- aft_predictions_individual(X, pred_times, m@theta,
                                       distribution = "log-normal")
  expect_equal(nrow(preds), n)
  expect_equal(ncol(preds), 3)
  # All predictions should be valid probabilities
  expect_true(all(unlist(preds) >= 0 & unlist(preds) <= 1))
})

test_that("aft_predictions_individual supports different distributions", {
  X <- cbind(1, 0)
  theta <- c(2, 0.3, log(1.5))
  times <- c(1, 5)

  # Weibull
  r1 <- aft_predictions_individual(X, times, theta,
                                    distribution = "weibull")
  expect_equal(nrow(r1), 1)

  # Log-logistic
  r2 <- aft_predictions_individual(X, times, theta,
                                    distribution = "log-logistic")
  expect_equal(nrow(r2), 1)

  # Log-normal
  r3 <- aft_predictions_individual(X, times, theta,
                                    distribution = "log-normal")
  expect_equal(nrow(r3), 1)
})

# aft_predictions_function tests (bd-9ejm.1) --------------------------------
#
# Function-level AFT predictions with delta-method point-wise confidence
# intervals, the R port of Python's utilities.aft_predictions_function.
#
# Expected contract (mirrors survival_predictions, which returns a data
# frame): aft_predictions_function(X, times, theta, covariance,
# distribution, measure = "survival", alpha = 0.05) returns a data frame
# with one row per time point and columns `time`, `predicted`, `variance`,
# `lower`, `upper`. `X` must be a single covariate pattern (one row); more
# than one row is an error whose message names the row/pattern constraint.
# `alpha` outside (0, 1) is an error.
#
# Fixture provenance: fixtures/aft_predictions_function.json, produced by
# fixtures/generate_aft_predictions_fixtures.py. That script fits AFT models
# with ee_aft on the deterministic Collett (2015) breast cancer data (no
# RNG) for the weibull, log-normal, log-logistic, and exponential
# distributions, then records aft_predictions_function output (predicted
# measure, delta-method variance, and Wald confidence limits) across five
# measures, two covariate patterns, and two alpha levels.
#
# Tolerances: point predictions are analytic and are pinned at 1e-6 (they
# agree with Python to ~1e-15). The delta-method variance is compared on
# the SE (square-root) scale, and the confidence limits directly, both at
# 1e-5. The looser 1e-5 is used because Python computes the delta-method
# variance with exact forward-mode automatic differentiation
# (delta_method(deriv_method = "exact")), whereas the R post-processing
# uses the delta_method default of central-difference numerical
# differentiation ("capprox"). The two agree to under 1e-6 relative on this
# fixture (~4.7e-7 relative on the SE scale); 1e-5 leaves margin for that
# documented method difference. The SE scale is used for the variance so the
# comparison stays in all.equal's relative mode even for the small-variance
# measures (see the inline note at the assertion).

test_that("aft_predictions_function matches Python across distributions and measures", {
  ref <- load_fixture("aft_predictions_function")

  for (i in seq_len(nrow(ref$scenarios))) {
    s <- ref$scenarios[i, ]
    dist <- s$distribution
    measure <- s$measure
    theta <- ref$models[[dist]]$theta
    covariance <- as.matrix(ref$models[[dist]]$covariance)
    X <- ref$scenarios$X[[i]]
    times <- ref$scenarios$times[[i]]
    alpha <- s$alpha

    result <- aft_predictions_function(
      X = X, times = times, theta = theta, covariance = covariance,
      distribution = dist, measure = measure, alpha = alpha
    )

    lbl <- paste0(dist, "/", measure, "/", s$pattern, "/alpha=", alpha)

    expect_equal(
      result$time, times,
      tolerance = 1e-6, label = paste0(lbl, ": time")
    )
    expect_equal(
      result$predicted, ref$scenarios$predicted[[i]],
      tolerance = 1e-6, label = paste0(lbl, ": predicted")
    )
    # Compare the variance on the SE (square-root) scale. testthat's
    # expect_equal follows all.equal, which switches to absolute comparison
    # when mean(abs(target)) <= tolerance. Several measures have variances
    # small enough (density down to ~1e-8, hazard to ~8e-7) that a direct
    # comparison against the 1e-5 tolerance would fall into absolute mode,
    # where a variance of 0 would pass vacuously. The SE vectors have mean
    # magnitude of at least ~1e-4, so every comparison stays relative.
    expect_equal(
      sqrt(result$variance), sqrt(ref$scenarios$variance[[i]]),
      tolerance = 1e-5, label = paste0(lbl, ": se")
    )
    expect_equal(
      result$lower, ref$scenarios$lower[[i]],
      tolerance = 1e-5, label = paste0(lbl, ": lower")
    )
    expect_equal(
      result$upper, ref$scenarios$upper[[i]],
      tolerance = 1e-5, label = paste0(lbl, ": upper")
    )
  }
})

test_that("aft_predictions_function with deriv_method = 'exact' matches Python", {
  # Python's aft_predictions_function computes the delta-method variance with
  # exact forward-mode autodiff. The R port forwards deriv_method to
  # delta_method; deriv_method = "exact" selects that same autodiff path, so it
  # reproduces the Python fixture to machine precision rather than to the looser
  # numerical-vs-exact tolerance used by the default capprox path above. The
  # variance is compared on the SE scale (see the section header for why a
  # direct comparison would fall into all.equal's vacuous absolute mode for the
  # small-variance measures); SE and confidence limits are pinned at 1e-8, which
  # both derivative implementations clear with margin (rounding only).
  ref <- load_fixture("aft_predictions_function")

  for (i in seq_len(nrow(ref$scenarios))) {
    s <- ref$scenarios[i, ]
    dist <- s$distribution
    measure <- s$measure
    theta <- ref$models[[dist]]$theta
    covariance <- as.matrix(ref$models[[dist]]$covariance)
    X <- ref$scenarios$X[[i]]
    times <- ref$scenarios$times[[i]]
    alpha <- s$alpha

    result <- aft_predictions_function(
      X = X, times = times, theta = theta, covariance = covariance,
      distribution = dist, measure = measure, alpha = alpha,
      deriv_method = "exact"
    )

    lbl <- paste0(dist, "/", measure, "/", s$pattern, "/alpha=", alpha)

    expect_equal(
      result$predicted, ref$scenarios$predicted[[i]],
      tolerance = 1e-6, label = paste0(lbl, ": predicted")
    )
    expect_equal(
      sqrt(result$variance), sqrt(ref$scenarios$variance[[i]]),
      tolerance = 1e-8, label = paste0(lbl, ": se")
    )
    expect_equal(
      result$lower, ref$scenarios$lower[[i]],
      tolerance = 1e-8, label = paste0(lbl, ": lower")
    )
    expect_equal(
      result$upper, ref$scenarios$upper[[i]],
      tolerance = 1e-8, label = paste0(lbl, ": upper")
    )
    expect_true(
      all(result$variance > 0),
      label = paste0(lbl, ": all variances strictly positive")
    )
  }
})

test_that("aft_predictions_function returns the expected data frame structure", {
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)
  times <- c(50, 100, 150, 200)

  result <- aft_predictions_function(
    X = matrix(c(1, 0), nrow = 1), times = times, theta = theta,
    covariance = covariance, distribution = "weibull", measure = "survival"
  )

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(times))
  expect_named(result, c("time", "predicted", "variance", "lower", "upper"))
})

test_that("aft_predictions_function point-wise CI brackets the prediction", {
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["log-normal"]]$theta
  covariance <- as.matrix(ref$models[["log-normal"]]$covariance)

  result <- aft_predictions_function(
    X = matrix(c(1, 1), nrow = 1), times = c(50, 100, 150, 200),
    theta = theta, covariance = covariance, distribution = "log-normal",
    measure = "risk"
  )

  expect_true(all(result$lower <= result$predicted))
  expect_true(all(result$upper >= result$predicted))
  expect_true(all(result$variance >= 0))
})

test_that("aft_predictions_function widens the interval with a smaller alpha", {
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)
  args <- list(
    X = matrix(c(1, 0), nrow = 1), times = c(50, 100, 150, 200),
    theta = theta, covariance = covariance, distribution = "weibull",
    measure = "risk"
  )

  ci_95 <- do.call(aft_predictions_function, c(args, list(alpha = 0.05)))
  ci_99 <- do.call(aft_predictions_function, c(args, list(alpha = 0.01)))

  expect_true(all((ci_99$upper - ci_99$lower) > (ci_95$upper - ci_95$lower)))
})

test_that("aft_predictions_function errors on multiple covariate patterns", {
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)

  # Contract for bd-9ejm.2: reject more than one covariate pattern, and the
  # message must name the row/pattern constraint (so this does not pass
  # vacuously on an object-not-found error before implementation).
  expect_error(
    aft_predictions_function(
      X = matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE),
      times = c(50, 100), theta = theta, covariance = covariance,
      distribution = "weibull"
    ),
    "row|pattern|single"
  )
})

test_that("aft_predictions_function errors on invalid alpha", {
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)

  expect_error(
    aft_predictions_function(
      X = matrix(c(1, 0), nrow = 1), times = 50, theta = theta,
      covariance = covariance, distribution = "weibull", alpha = 0
    ),
    "alpha"
  )
  expect_error(
    aft_predictions_function(
      X = matrix(c(1, 0), nrow = 1), times = 50, theta = theta,
      covariance = covariance, distribution = "weibull", alpha = 1
    ),
    "alpha"
  )
})

test_that("aft_predictions_function handles a length-one times vector", {
  # R improvement over Python, which raises on a scalar/length-one `times`.
  # A single-time call must equal the matching row of a multi-time call.
  ref <- load_fixture("aft_predictions_function")
  theta <- ref$models[["weibull"]]$theta
  covariance <- as.matrix(ref$models[["weibull"]]$covariance)
  args <- list(
    X = matrix(c(1, 0), nrow = 1), theta = theta, covariance = covariance,
    distribution = "weibull", measure = "risk"
  )

  many <- do.call(
    aft_predictions_function, c(args, list(times = c(50, 100, 150, 200)))
  )
  one <- do.call(aft_predictions_function, c(args, list(times = 100)))

  expect_equal(nrow(one), 1)
  expect_named(one, c("time", "predicted", "variance", "lower", "upper"))
  expect_equal(one[1, ], many[many$time == 100, ], ignore_attr = TRUE)
})

# plogit_predict tests -------------------------------------------------------

test_that("plogit_predict returns correct shape with default S", {
  set.seed(42)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  # Generate discrete survival data
  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  # Fit the model
  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

 result <- plogit_predict(m@theta, time = y, event = delta, X = X)

  # Should be K-by-n matrix
  expect_true(is.matrix(result))
  expect_equal(nrow(result), k)
  expect_equal(ncol(result), n)
})

test_that("plogit_predict survival is between 0 and 1", {
  set.seed(123)
  n <- 100
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  # Generate discrete survival data
  max_time <- 6
  hazard_prob <- inverse_logit(-3 + 0.5 * x)
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  for (i in seq_len(n)) {
    for (k_t in seq_len(max_time)) {
      if (runif(1) < hazard_prob[i]) {
        y[i] <- k_t
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  result <- plogit_predict(m@theta, time = y, event = delta, X = X)

  # Survival should be in [0, 1]
  expect_true(all(result >= 0 & result <= 1))
})

test_that("plogit_predict survival decreases over time", {
  set.seed(456)
  n <- 100
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 8
  hazard_prob <- inverse_logit(-3 + 0.3 * x)
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  for (i in seq_len(n)) {
    for (k_t in seq_len(max_time)) {
      if (runif(1) < hazard_prob[i]) {
        y[i] <- k_t
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  result <- plogit_predict(m@theta, time = y, event = delta, X = X,
                            measure = "survival")

  # For each individual, survival should be non-increasing
  for (j in seq_len(ncol(result))) {
    expect_true(all(diff(result[, j]) <= 1e-10))
  }
})

test_that("plogit_predict risk equals 1 - survival", {
  set.seed(789)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  surv <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "survival")
  risk <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "risk")

  expect_equal(risk, 1 - surv)
})

test_that("plogit_predict cumulative_hazard equals -log(survival)", {
  set.seed(101)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  surv <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "survival")
  chaz <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "cumulative_hazard")

  expect_equal(chaz, -log(surv))
})

test_that("plogit_predict with times_to_predict subsets correctly", {
  set.seed(202)
  n <- 100
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 8
  hazard_prob <- inverse_logit(-3 + 0.3 * x)
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  for (i in seq_len(n)) {
    for (k_t in seq_len(max_time)) {
      if (runif(1) < hazard_prob[i]) {
        y[i] <- k_t
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  # Predict at specific time points
  pred_times <- c(1, 3, 5)
  result <- plogit_predict(m@theta, time = y, event = delta, X = X,
                            times_to_predict = pred_times)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), n)
})

test_that("plogit_predict at time 0 returns baseline", {
  set.seed(303)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  # At time 0, survival should be 1 for all individuals
  result <- plogit_predict(m@theta, time = y, event = delta, X = X,
                            times_to_predict = c(0),
                            measure = "survival")
  expect_equal(as.numeric(result), rep(1, n))
})

test_that("plogit_predict errors beyond max observed time", {
  set.seed(404)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  expect_error(
    plogit_predict(m@theta, time = y, event = delta, X = X,
                    times_to_predict = c(100)),
    "maximum observed time"
  )
})

test_that("plogit_predict works with custom S matrix", {
  set.seed(606)
  n <- 200
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 8
  y <- rep(NA_real_, n)
  delta <- rep(0, n)
  for (i in seq_len(n)) {
    for (k_t in seq_len(max_time)) {
      h <- inverse_logit(-4 + 0.3 * x[i] + 0.1 * k_t)
      if (runif(1) < h) {
        y[i] <- k_t
        delta[i] <- 1
        break
      }
    }
    if (is.na(y[i])) {
      y[i] <- max_time
      delta[i] <- 0
    }
  }

  # Gompertz-like S: intercept + linear time
  t_steps <- seq_len(max_time)
  S_mat <- cbind(1, t_steps)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta, S = S_mat)
  }
  inits <- c(0, -4, 0.1)
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  result <- plogit_predict(m@theta, time = y, event = delta, X = X,
                            S = S_mat, measure = "survival")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), max_time)
  expect_equal(ncol(result), n)
  # Survival should be in [0, 1]
  expect_true(all(result >= 0 & result <= 1))
})

test_that("plogit_predict density equals hazard * survival", {
  set.seed(505)
  n <- 50
  x <- rbinom(n, 1, 0.5)
  X <- cbind(x)

  max_time <- 5
  y <- sample(1:max_time, n, replace = TRUE)
  delta <- rbinom(n, 1, 0.6)

  unique_event_times <- sort(unique(y[delta == 1]))
  k <- length(unique_event_times)

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = y, event = delta)
  }
  inits <- c(0, -3, rep(0, k - 1))
  m <- MEstimator(stacked_equations = psi, init = inits)
  m <- estimate(m, solver = "nleqslv")

  surv <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "survival")
  haz <- plogit_predict(m@theta, time = y, event = delta, X = X,
                         measure = "hazard")
  dens <- plogit_predict(m@theta, time = y, event = delta, X = X,
                          measure = "density")

  expect_equal(dens, haz * surv)
})

# ---- consistent time/event argument names (bd-9ejm.4 contract) --------------
#
# See the header in test-ee-survival.R for the full naming contract. These
# tests enforce that plogit_predict uses
#
#   plogit_predict(theta, time, event, X, S, times_to_predict, measure,
#                  unique_times)
#
# so that `time` and `event` match ee_plogit and the other survival functions.
# The prediction-grid argument `times_to_predict` is a separate concept and is
# left unchanged. The previous signature named these arguments `t` and `delta`,
# and the rejection test below asserts that those names no longer exist:
#
#   plogit_predict(theta, t, delta, X, S, times_to_predict, measure,
#                  unique_times)
#
# The plogit_predict fixture pins the ee_plogit -> plogit_predict pipeline end
# to end with the Python semantics (t = observed time, delta = event). It was
# generated by generate_plogit_predict_fixture() in
# fixtures/generate_all_remaining_fixtures.py, which fits ee_plogit on the
# breast cancer data with a Gompertz-like time specification and records the
# predicted risk at several interior times.

test_that("plogit_predict accepts time/event and matches Python fixture", {
  ref <- load_fixture("plogit_predict")
  X <- ref$X
  S <- ref$S

  result <- plogit_predict(
    ref$theta,
    time = ref$t, event = ref$delta, X = X, S = S,
    times_to_predict = ref$times_to_predict, measure = ref$measure
  )

  expect_equal(result, as.matrix(ref$predictions), tolerance = 1e-6)
})

test_that("plogit_predict rejects the old t/delta argument names", {
  ref <- load_fixture("plogit_predict")

  # `t` was the observed time and `delta` the event indicator; after the
  # rename neither name exists.
  expect_error(
    plogit_predict(ref$theta, t = ref$t, delta = ref$delta, X = ref$X,
                   S = ref$S, times_to_predict = ref$times_to_predict,
                   measure = ref$measure)
  )
})

test_that("ee_plogit + plogit_predict with time/event agree end to end", {
  # Fit with ee_plogit using time=/event=, then predict with plogit_predict
  # using the SAME named data. This is the cross-function consistency check:
  # both functions must share one meaning of `time` and `event`, and the
  # result must reproduce the Python reference.
  ref <- load_fixture("plogit_predict")
  X <- ref$X
  S <- ref$S

  psi <- function(theta) {
    ee_plogit(theta, X = X, time = ref$t, event = ref$delta, S = S)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m)

  expect_equal(unname(m@theta), ref$theta, tolerance = 1e-4)

  result <- plogit_predict(
    m@theta,
    time = ref$t, event = ref$delta, X = X, S = S,
    times_to_predict = ref$times_to_predict, measure = ref$measure
  )

  expect_equal(result, as.matrix(ref$predictions), tolerance = 1e-4)
})

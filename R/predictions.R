#' Generate predicted values from a regression model
#'
#' Computes predicted outcomes, their variance, and Wald-type confidence
#' intervals from estimated regression coefficients and their covariance
#' matrix. This is a post-processing utility meant to be used after
#' [MEstimator()] has been fitted.
#'
#' No transformations are applied. For logistic models this returns log-odds
#' (not probabilities). Apply [stats::plogis()] for the probability scale, or
#' [inverse_logit()] if the values feed a transform passed to [delta_method()]
#' with `deriv_method = "exact"`.
#'
#' @param X Numeric n-by-p design matrix of covariate values for prediction.
#' @param theta Numeric vector of p estimated coefficients (from
#'   `m@theta`).
#' @param covariance Numeric p-by-p covariance matrix (from `m@variance`).
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#' @param alpha Numeric significance level for confidence intervals.
#'   Default `0.05` (95% CIs).
#'
#' @returns A data frame with n rows and columns: `predicted`, `variance`,
#'   `lower`, `upper`.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.5))
#' dat$y <- 1 + 0.5 * dat$x + 2 * dat$z + rnorm(n)
#'
#' m <- m_estimate(y ~ x + z, data = dat, .ee = ee_regression, model = "linear")
#'
#' # Predict along a small grid of x, holding z at 1. The columns of the grid
#' # must match the order of the coefficients, intercept first.
#' X_new <- cbind(1, x = c(-1, 0, 1), z = 1)
#' regression_predictions(X_new, theta = coef(m), covariance = vcov(m))
#'
#' @export
regression_predictions <- function(
  X,
  theta,
  covariance,
  offset = NULL,
  alpha = 0.05
) {
  check_prediction_alpha(alpha)

  X <- as.matrix(X)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  # Predicted Y = X %*% theta + offset
  yhat <- as.numeric(X %*% theta)
  if (!is.null(offset)) {
    yhat <- yhat + as.numeric(offset)
  }

  # Predicted Y variance: diag(X %*% Sigma %*% t(X))
  # Efficient: rowSums((X %*% Sigma) * X)
  yhat_var <- rowSums((X %*% covariance) * X)

  # Confidence intervals. The variance is a quadratic form in a covariance
  # matrix rather than a numerical derivative, so it is not floored at zero the
  # way the two delta-method helpers below floor theirs.
  band <- wald_band(yhat, sqrt(yhat_var), alpha)

  data.frame(
    predicted = yhat,
    variance = yhat_var,
    lower = band$lower,
    upper = band$upper
  )
}

#' Generate predicted survival measures from a parametric survival model
#'
#' Computes predicted survival analysis measures and point-wise confidence
#' intervals using the delta method. Meant to be used after fitting
#' [ee_survival_model()] with [MEstimator()].
#'
#' @param times Numeric vector of time points for prediction.
#' @param theta Numeric vector of estimated parameters from
#'   `ee_survival_model`.
#' @param covariance Numeric covariance matrix from `m@variance`.
#' @param distribution Character string matching the distribution used in
#'   `ee_survival_model`: `"exponential"`, `"weibull"`, or `"gompertz"`. Any
#'   other value is an error.
#' @param measure Character string: `"survival"`, `"risk"`, `"density"`,
#'   `"hazard"`, or `"cumulative_hazard"`. Default `"survival"`.
#' @param alpha Numeric significance level. Default `0.05`.
#' @param deriv_method Character string for the derivative method used to build
#'   the delta-method Jacobian. One of `"capprox"` (central difference),
#'   `"fapprox"` (forward difference), `"bapprox"` (backward difference), or
#'   `"exact"` (forward-mode automatic differentiation). Default `"capprox"`.
#'   Python Delicatessen uses exact differentiation internally; pass
#'   `deriv_method = "exact"` to reproduce it with exact derivatives and no
#'   step-size tuning. See [delta_method()].
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. Default `1e-9`. The step is absolute and is
#'   floored at the floating-point resolution of each estimate, so a large
#'   parameter magnitude cannot silently reduce it to nothing; see
#'   [approx_differentiation()].
#'
#' @returns A data frame with columns: `time`, `predicted`, `variance`,
#'   `lower`, `upper`.
#'
#' @section Gompertz distribution:
#' The Gompertz survival and hazard follow the [ee_survival_model()]
#' parameterization,
#'
#' \deqn{S(t) = \exp\left(-\frac{\lambda}{\gamma}\left(e^{\gamma t} - 1\right)
#'   \right), \qquad h(t) = \lambda e^{\gamma t}.}
#'
#' This is a deliberate divergence from Python Delicatessen, whose
#' `survival_predictions` branches only on the exponential distribution and
#' otherwise applies the Weibull formulas, so a `"gompertz"` request there
#' silently returns Weibull values.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' # Weibull survival times for 45 women with breast cancer, with no covariates.
#' # The default rootSolve solver does not converge here, so nleqslv is used.
#' psi <- function(theta) {
#'   ee_survival_model(
#'     theta,
#'     time = breast_cancer$times,
#'     event = breast_cancer$delta,
#'     distribution = "weibull"
#'   )
#' }
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(0.1, 0.1),
#'   solver = "nleqslv"
#' )
#'
#' # The survival function at three follow-up times, with delta-method
#' # confidence intervals. Observed times run from 5 to 225 days.
#' survival_predictions(
#'   times = c(50, 100, 150),
#'   theta = m@theta,
#'   covariance = m@variance,
#'   distribution = "weibull",
#'   measure = "survival"
#' )
#'
#' @export
survival_predictions <- function(
  times,
  theta,
  covariance,
  distribution,
  measure = "survival",
  alpha = 0.05,
  deriv_method = "capprox",
  dx = 1e-9
) {
  check_prediction_alpha(alpha)

  distribution <- tolower(distribution)
  times <- as.numeric(times)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  # Validate the distribution up front so an unsupported choice aborts before
  # delta_method runs, rather than silently applying the Weibull formulas.
  supported <- c("exponential", "weibull", "gompertz")
  if (!distribution %in% supported) {
    cli::cli_abort(
      c(
        "The distribution {.val {distribution}} is not supported.",
        "i" = "Use one of: {.val exponential}, {.val weibull}, {.val gompertz}."
      )
    )
  }

  # Function to compute survival metric for a single time point. The scalar
  # arithmetic keeps derivatives attached under exact autodiff, where `th` is a
  # tangent-carrying pair vector.
  predict_at_time <- function(t_val, th) {
    lambd <- th[1]
    if (distribution == "exponential") {
      gamma <- 1
    } else {
      gamma <- th[2]
    }
    if (distribution == "gompertz") {
      # Gompertz parameterization matching ee_survival_model.
      surv <- exp(-(lambd / gamma) * (exp(gamma * t_val) - 1))
      haz <- lambd * exp(gamma * t_val)
    } else {
      surv <- exp(-lambd * t_val^gamma)
      haz <- lambd * gamma * t_val^(gamma - 1)
    }
    convert_survival_measures(surv, haz, measure)
  }

  # Point estimates for each time
  est <- vapply(times, function(t) predict_at_time(t, theta), numeric(1))

  # Delta method for variance at each time. The transform concatenates the
  # per-time measures with c() rather than coercing with vapply(numeric(1)):
  # under exact autodiff the per-time values are tangent-carrying pair objects,
  # not plain doubles, so a vapply(numeric(1)) template cannot hold them.
  g_func <- function(th) {
    do.call(c, lapply(times, function(t) predict_at_time(t, th)))
  }
  cov_m <- delta_method(
    theta,
    transform = g_func,
    covariance = covariance,
    deriv_method = deriv_method,
    dx = dx
  )
  var_m <- diag(cov_m)

  # Confidence intervals. pmax(var_m, 0) guards the finite-difference paths,
  # where step-size cancellation can drive a near-zero variance slightly
  # negative; it is a no-op under exact autodiff, whose delta-method variance
  # is nonnegative by construction. Python does not clamp.
  band <- wald_band(est, sqrt(pmax(var_m, 0)), alpha)

  data.frame(
    time = times,
    predicted = est,
    variance = var_m,
    lower = band$lower,
    upper = band$upper
  )
}

#' Predicted survival measures from an AFT model
#'
#' Computes individual-level predicted survival measures from an
#' accelerated failure time model at specified time points. Meant to be
#' used after fitting [ee_aft()] with [MEstimator()].
#'
#' @param X Numeric n-by-b design matrix of covariate values.
#' @param times Numeric vector of time points for prediction.
#' @param theta Numeric vector of estimated parameters from `ee_aft`.
#' @param distribution Character string matching the distribution used in
#'   `ee_aft`.
#' @param measure Character string: `"survival"`, `"risk"`, `"density"`,
#'   `"hazard"`, or `"cumulative_hazard"`. Default `"survival"`.
#'
#' @returns A data frame with n rows and one column per time point.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' # Weibull AFT fit, then individual-level survival for the first four people
#' set.seed(1)
#' n <- 200
#' x <- rbinom(n, 1, 0.5)
#' Xd <- cbind(1, x)
#' eps <- log(-log(runif(n)))
#' t_event <- exp(2 + 0.5 * x + 0.8 * eps)
#' t_censor <- rexp(n, rate = 0.02)
#' t_obs <- pmin(t_event, t_censor)
#' delta <- as.numeric(t_event <= t_censor)
#'
#' psi <- function(theta) {
#'   ee_aft(theta, X = Xd, time = t_obs, event = delta, distribution = "weibull")
#' }
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(mean(log(t_obs)), 0, 0),
#'   solver = "nleqslv"
#' )
#'
#' # Rows are individuals and columns are the requested times. The first two
#' # people share a covariate pattern, as do the third and fourth, so their
#' # predictions agree.
#' aft_predictions_individual(
#'   X = Xd[1:4, ], times = c(5, 10, 20),
#'   theta = m@theta, distribution = "weibull", measure = "survival"
#' )
#'
#' @export
aft_predictions_individual <- function(
  X,
  times,
  theta,
  distribution,
  measure = "survival"
) {
  X <- as.matrix(X)
  theta <- as.numeric(theta)
  distribution <- tolower(distribution)
  n <- nrow(X)
  beta_dim <- ncol(X)
  beta <- theta[seq_len(beta_dim)]

  # Extract sigma: exp(-theta_last) for non-exponential
  if (distribution == "exponential") {
    sigma <- 1
  } else {
    sigma <- exp(-theta[length(theta)])
  }

  # Linear predictor, one entry per person. The matrix product is what makes
  # this the vectorized half of the pair: aft_measure_at_time() is elementwise
  # throughout, so a whole column of people is predicted in one call.
  xbeta <- as.numeric(X %*% beta)

  # Compute survival metric at each time point
  results <- matrix(NA_real_, nrow = n, ncol = length(times))

  for (j in seq_along(times)) {
    results[, j] <- aft_measure_at_time(
      times[j],
      xbeta,
      sigma,
      distribution,
      measure
    )
  }

  # Return as data frame with time-named columns
  colnames(results) <- paste0("t_", times)
  as.data.frame(results)
}

#' Function-level predicted survival measures from an AFT model
#'
#' Computes predicted survival analysis measures and point-wise confidence
#' intervals from an accelerated failure time model for a single covariate
#' pattern across a set of time points. The point estimates mirror
#' [aft_predictions_individual()]; the variance is obtained with the delta
#' method (see [delta_method()]) and Wald-type intervals are formed on the
#' resulting standard errors. Meant to be used after fitting [ee_aft()] with
#' [MEstimator()], typically to draw a measure and its confidence band over
#' time.
#'
#' The survival and hazard for a covariate pattern \eqn{X} are
#'
#' \deqn{S(t) = S_{\epsilon}\left( \frac{\log(t) - X \beta^T}{\sigma} \right)}
#' \deqn{h(t) = (\sigma t)^{-1} h_{\epsilon}\left(
#'   \frac{\log(t) - X \beta^T}{\sigma} \right)}
#'
#' where \eqn{S_{\epsilon}} and \eqn{h_{\epsilon}} are the error survival and
#' hazard functions for the chosen `distribution`. The requested `measure` is
#' derived from these through [convert_survival_measures()].
#'
#' @param X Numeric 1-by-b design matrix giving a single covariate pattern.
#'   More than one row is an error, since each pattern has its own variance.
#' @param times Numeric vector of time points for prediction.
#' @param theta Numeric vector of estimated parameters from `ee_aft`.
#' @param covariance Numeric covariance matrix from `m@variance`.
#' @param distribution Character string matching the distribution used in
#'   `ee_aft`.
#' @param measure Character string: `"survival"`, `"risk"`, `"density"`,
#'   `"hazard"`, or `"cumulative_hazard"`. Default `"survival"`.
#' @param alpha Numeric significance level. Default `0.05` (95% CIs).
#' @param deriv_method Character string for the derivative method used to build
#'   the delta-method Jacobian. One of `"capprox"` (central difference),
#'   `"fapprox"` (forward difference), `"bapprox"` (backward difference), or
#'   `"exact"` (forward-mode automatic differentiation). Default `"capprox"`.
#'   Python Delicatessen uses exact differentiation internally; pass
#'   `deriv_method = "exact"` to reproduce it with exact derivatives and no
#'   step-size tuning. See [delta_method()].
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. Default `1e-9`. The step is absolute and is
#'   floored at the floating-point resolution of each estimate, so a large
#'   parameter magnitude cannot silently reduce it to nothing; see
#'   [approx_differentiation()].
#'
#' @returns A data frame with one row per time point and columns: `time`,
#'   `predicted`, `variance`, `lower`, `upper`.
#'
#' @section Length-one times:
#' A single time point is supported and returns a one-row data frame equal to
#' the corresponding row of a multi-time call. This is a deliberate
#' improvement over Python Delicatessen, whose `aft_predictions_function`
#' raises on a scalar or length-one `times` because it takes the diagonal of a
#' scalar delta-method covariance.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' # Weibull AFT fit, then a survival curve for one covariate pattern
#' set.seed(1)
#' n <- 200
#' x <- rbinom(n, 1, 0.5)
#' Xd <- cbind(1, x)
#' eps <- log(-log(runif(n)))
#' t_event <- exp(2 + 0.5 * x + 0.8 * eps)
#' t_censor <- rexp(n, rate = 0.02)
#' t_obs <- pmin(t_event, t_censor)
#' delta <- as.numeric(t_event <= t_censor)
#'
#' psi <- function(theta) {
#'   ee_aft(theta, X = Xd, time = t_obs, event = delta, distribution = "weibull")
#' }
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(mean(log(t_obs)), 0, 0),
#'   solver = "nleqslv"
#' )
#'
#' aft_predictions_function(
#'   X = matrix(c(1, 1), nrow = 1), times = c(5, 10, 20),
#'   theta = m@theta, covariance = m@variance,
#'   distribution = "weibull", measure = "risk"
#' )
#'
#' @export
aft_predictions_function <- function(
  X,
  times,
  theta,
  covariance,
  distribution,
  measure = "survival",
  alpha = 0.05,
  deriv_method = "capprox",
  dx = 1e-9
) {
  check_prediction_alpha(alpha)

  X <- as.matrix(X)

  # A single covariate pattern is required: each pattern has its own variance,
  # so more than one row is rejected rather than tracked jointly.
  if (nrow(X) > 1) {
    cli::cli_abort(
      c(
        "{.arg X} must be a single covariate pattern (one row).",
        "x" = "{nrow(X)} rows were provided."
      )
    )
  }

  distribution <- tolower(distribution)
  times <- as.numeric(times)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  x_vec <- as.numeric(X[1, ])
  beta_dim <- length(x_vec)
  n_theta <- length(theta)

  # Measure for the single covariate pattern at one time point. The linear
  # predictor is accumulated one term at a time rather than through a matrix
  # product, because that is what lets exact differentiation reach here: under
  # exact autodiff `th` is a tangent-carrying pair vector, and indexing plus
  # scalar arithmetic keep the derivatives attached, where a matrix product
  # would need as.numeric() and a tangent-carrying value cannot become a plain
  # double. The same code returns plain doubles when `th` is numeric, so the
  # finite-difference paths reuse it.
  predict_at_time <- function(t_val, th) {
    xbeta <- 0
    for (k in seq_len(beta_dim)) {
      xbeta <- xbeta + x_vec[k] * th[k]
    }
    if (distribution == "exponential") {
      sigma <- 1
    } else {
      sigma <- exp(-th[n_theta])
    }
    aft_measure_at_time(t_val, xbeta, sigma, distribution, measure)
  }

  # Measure across times for the single covariate pattern. The transform
  # concatenates the per-time measures with c() rather than coercing with
  # vapply(numeric(1)): under exact autodiff the per-time values are tangent-
  # carrying pair objects, not plain doubles, so a vapply(numeric(1)) template
  # cannot hold them.
  g_func <- function(th) {
    do.call(c, lapply(times, function(t) predict_at_time(t, th)))
  }

  # Point estimates for each time (numeric theta, no tangents present)
  est <- vapply(times, function(t) predict_at_time(t, theta), numeric(1))

  # Delta-method covariance of the measure across times, then its diagonal
  cov_m <- delta_method(
    theta,
    transform = g_func,
    covariance = covariance,
    deriv_method = deriv_method,
    dx = dx
  )
  var_m <- diag(as.matrix(cov_m))

  # Wald-type point-wise confidence intervals. pmax(var_m, 0) guards the finite-
  # difference paths, where step-size cancellation can drive a near-zero
  # variance slightly negative; it is a no-op under exact autodiff. Python does
  # not clamp.
  band <- wald_band(est, sqrt(pmax(var_m, 0)), alpha)

  data.frame(
    time = times,
    predicted = est,
    variance = var_m,
    lower = band$lower,
    upper = band$upper
  )
}

#' Predicted survival measures from a pooled logistic regression model
#'
#' Computes predicted survival analysis measures from a pooled logistic
#' regression model at specified time points. Meant to be used after
#' fitting [ee_plogit()] with [MEstimator()].
#'
#' @param theta Numeric vector of estimated parameters from `ee_plogit`.
#' @param time Numeric vector of n observed (possibly censored) times (same
#'   as in `ee_plogit`).
#' @param event Numeric vector of n event indicators (same as in
#'   `ee_plogit`).
#' @param X Numeric n-by-b design matrix for covariates.
#' @param S Optional time design matrix. Default `NULL` uses disjoint
#'   indicators.
#' @param times_to_predict Optional numeric vector of specific times to
#'   predict at. Default `NULL` returns all time steps.
#' @param measure Character string: `"survival"`, `"risk"`, `"density"`,
#'   `"hazard"`, or `"cumulative_hazard"`. Default `"survival"`.
#' @param unique_times Optional numeric vector of unique event times.
#'   Default `NULL`.
#'
#' @returns A K-by-n matrix (or selected subset of rows if
#'   `times_to_predict` is specified).
#'
#' @examples
#' # Bladder tumor recurrence, fit with disjoint indicators for time
#' d <- collett_bladder
#' d$novel <- d$treat - 1
#' W <- as.matrix(d[, c("novel", "init", "size")])
#' k <- length(unique(d$time[d$delta == 1]))
#'
#' psi <- function(theta) {
#'   ee_plogit(theta, X = W, time = d$time, event = d$delta)
#' }
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(rep(0, ncol(W)), -4, rep(0, k - 1))
#' )
#'
#' # Rows are the requested times and columns are individuals, so this shows
#' # disease-free survival at 12, 24, and 59 months for the first five people.
#' plogit_predict(
#'   m@theta,
#'   time = d$time, event = d$delta, X = W,
#'   times_to_predict = c(12, 24, 59), measure = "survival"
#' )[, 1:5]
#'
#' @export
plogit_predict <- function(
  theta,
  time,
  event,
  X,
  S = NULL,
  times_to_predict = NULL,
  measure = "survival",
  unique_times = NULL
) {
  time <- as.numeric(time)
  event <- as.numeric(event)
  X <- as.matrix(X)
  n <- nrow(X)
  xp <- ncol(X)

  # Split theta
  beta_x <- theta[seq_len(xp)]
  beta_s <- theta[(xp + 1):length(theta)]

  # Build time design matrix
  if (is.null(S)) {
    if (is.null(unique_times)) {
      event_times <- time[event == 1]
      unique_times <- sort(unique(event_times))
    } else {
      unique_times <- as.numeric(unique_times)
    }
    n_time_steps <- length(unique_times)
    # Disjoint indicator matrix with first column as intercept
    time_design_matrix <- diag(n_time_steps)
    time_design_matrix[, 1] <- 1
  } else {
    S <- as.matrix(S)
    time_design_matrix <- S
    unique_times <- seq_len(max(time))
    n_time_steps <- length(unique_times)
  }

  # Log-odds contributions
  log_odds_w <- as.numeric(X %*% beta_x) # n-vector
  log_odds_t <- as.numeric(time_design_matrix %*% beta_s) # K-vector

  # Predicted event probability at each time: K-by-n matrix
  log_odds_w_matrix <- matrix(
    rep(log_odds_w, each = n_time_steps),
    nrow = n_time_steps,
    ncol = n
  )
  y_pred <- inverse_logit(log_odds_w_matrix + log_odds_t)

  # Cumulative product of (1 - p) gives survival: S(t) = prod(1 - h(t_k)).
  # apply() drops to a bare vector when each column's cumprod has length one
  # (a single unique event time), so reshape back to K-by-n: a no-op for
  # K > 1 that restores the matrix contract for K = 1.
  survival_pred <- matrix(
    apply(1 - y_pred, 2, cumprod),
    nrow = n_time_steps,
    ncol = n
  )

  # Convert to requested measure
  prediction_matrix <- convert_survival_measures(
    survival_pred,
    hazard = y_pred,
    measure = measure
  )

  # Return predictions
  if (is.null(times_to_predict)) {
    return(prediction_matrix)
  }

  # Subset to requested time points
  prediction_t0 <- convert_survival_measures(1, 0, measure = measure)
  results <- matrix(NA_real_, nrow = length(times_to_predict), ncol = n)

  for (i in seq_along(times_to_predict)) {
    time_val <- times_to_predict[i]

    if (time_val == 0 || time_val < unique_times[1]) {
      # Before first event time
      results[i, ] <- rep(prediction_t0, n)
    } else if (time_val > max(time)) {
      cli::cli_abort("Cannot predict beyond the maximum observed time.")
    } else {
      # Find appropriate index in prediction matrix
      if (unique_times[n_time_steps] <= time_val) {
        idx <- n_time_steps
      } else {
        # Find the last unique time <= time_val
        idx <- max(which(unique_times <= time_val))
      }
      results[i, ] <- prediction_matrix[idx, ]
    }
  }

  results
}

# ---- Shared internals --------------------------------------------------------
# Three pieces of the functions above are the same piece, written once here so
# that a change to any of them reaches every function that uses it.
#
# `aft_measure_at_time()` is the one to read before changing anything. It is
# reached both with plain doubles and, under `deriv_method = "exact"`, with
# tangent-carrying pairs, and it is written in scalar arithmetic for the second
# case. Its two callers form the linear predictor differently and have to keep
# doing so: `aft_predictions_individual()` uses a matrix product, which predicts
# every person at once but which the exact path cannot use, because reading a
# matrix product back out needs `as.numeric()` and a tangent-carrying value
# cannot become a plain double. Keeping the linear predictor in the callers and
# the per-time measure here lets the individual predictions stay vectorized over
# people while the function-level predictions stay differentiable.

#' Validate a significance level for the prediction helpers
#'
#' These helpers replicate the Python Delicatessen API down to the message they
#' raise, which is why they do not share `check_alpha()`: that one also rejects
#' a non-numeric or non-scalar `alpha` and words its message differently. The
#' error is reported against the caller so that it names the prediction function
#' the user called rather than this helper.
#'
#' @param alpha The significance level supplied.
#' @returns Invisible `NULL`. Raises an error if `alpha` is not strictly
#'   between 0 and 1.
#' @noRd
check_prediction_alpha <- function(alpha) {
  if (alpha <= 0 || alpha >= 1) {
    cli::cli_abort(
      "{.arg alpha} must be between 0 and 1 (exclusive).",
      call = rlang::caller_env()
    )
  }
  invisible(NULL)
}

#' Wald confidence limits
#'
#' The symmetric limits `estimate` plus and minus the two-sided standard normal
#' critical value times `se`.
#'
#' The standard error is taken already formed rather than derived from a
#' variance, because the callers do not agree on how to reach it. The two
#' delta-method helpers floor the variance at zero first, since a
#' finite-difference derivative can drive a near-zero variance slightly
#' negative, while [regression_predictions()] computes a quadratic form and does
#' not floor it.
#'
#' @param estimate Numeric vector of point estimates.
#' @param se Numeric vector of standard errors.
#' @param alpha Numeric significance level.
#' @returns A list with numeric `lower` and `upper` elements, each the length of
#'   `estimate`.
#' @noRd
wald_band <- function(estimate, se, alpha) {
  z_alpha <- qnorm(1 - alpha / 2)
  list(
    lower = estimate - z_alpha * se,
    upper = estimate + z_alpha * se
  )
}

#' The predicted measure of an AFT model at one time
#'
#' Evaluates the survival and hazard of the error distribution at
#' \eqn{\epsilon = (\log t - X \beta) / \sigma} and converts them to the
#' requested measure. The hazard scaler \eqn{1 / (\sigma t)} carries the hazard
#' from the error scale to the time scale.
#'
#' `xbeta` and `sigma` arrive already formed, for the reason given in the
#' section comment above. Everything here is elementwise arithmetic plus
#' [standard_normal_cdf()] and [standard_normal_pdf()], all of which have exact
#' differentiation rules, so a tangent handed in through `xbeta` or `sigma`
#' comes back out in the result, and a vector `xbeta` of one entry per person
#' comes back as a vector of the same length. `pnorm()` and `dnorm()` would
#' serve the numeric callers equally but hand their argument to compiled code
#' without dispatching, so they cannot serve the exact path.
#'
#' @param t_val The time to predict at.
#' @param xbeta The linear predictor: a scalar, or one value per person.
#' @param sigma The scale parameter.
#' @param distribution The error distribution, already lowercased.
#' @param measure The measure to return.
#' @returns The requested measure, in the shape of `xbeta`.
#' @noRd
aft_measure_at_time <- function(t_val, xbeta, sigma, distribution, measure) {
  eps <- (log(t_val) - xbeta) / sigma
  hazard_scaler <- 1 / (sigma * t_val)

  if (distribution %in% c("exponential", "weibull")) {
    surv <- exp(-exp(eps))
    haz <- hazard_scaler * exp(eps)
  } else if (distribution %in% c("lognormal", "log-normal")) {
    surv <- 1 - standard_normal_cdf(eps)
    haz <- hazard_scaler * standard_normal_pdf(eps) / surv
  } else if (distribution %in% c("loglogistic", "log-logistic")) {
    surv <- 1 / (1 + exp(eps))
    haz <- hazard_scaler / (1 + exp(-eps))
  } else {
    cli::cli_abort(
      c(
        "The distribution {.val {distribution}} is not supported.",
        "i" = "Use one of: {.val exponential}, {.val weibull},
               {.val log-logistic}, {.val log-normal}."
      ),
      call = rlang::caller_env()
    )
  }
  convert_survival_measures(surv, haz, measure)
}

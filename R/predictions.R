#' Generate predicted values from a regression model
#'
#' Computes predicted outcomes, their variance, and Wald-type confidence
#' intervals from estimated regression coefficients and their covariance
#' matrix. This is a post-processing utility meant to be used after
#' [MEstimator()] has been fitted.
#'
#' No transformations are applied. For logistic models this returns
#' log-odds (not probabilities). Apply [inverse_logit()] to transform.
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
#' @export
regression_predictions <- function(X, theta, covariance,
                                   offset = NULL, alpha = 0.05) {
  # Validate alpha
  if (alpha <= 0 || alpha >= 1) {
    cli::cli_abort("{.arg alpha} must be between 0 and 1 (exclusive).")
  }

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

  # Confidence intervals
  yhat_se <- sqrt(yhat_var)
  z_alpha <- qnorm(1 - alpha / 2)
  lower_ci <- yhat - z_alpha * yhat_se
  upper_ci <- yhat + z_alpha * yhat_se

  data.frame(
    predicted = yhat,
    variance = yhat_var,
    lower = lower_ci,
    upper = upper_ci
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
#'   `ee_survival_model`.
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
#'   `deriv_method = "exact"`. Default `1e-9`.
#'
#' @returns A data frame with columns: `time`, `predicted`, `variance`,
#'   `lower`, `upper`.
#'
#' @export
survival_predictions <- function(times, theta, covariance, distribution,
                                  measure = "survival", alpha = 0.05,
                                  deriv_method = "capprox", dx = 1e-9) {
  if (alpha <= 0 || alpha >= 1) {
    cli::cli_abort("{.arg alpha} must be between 0 and 1 (exclusive).")
  }

  distribution <- tolower(distribution)
  times <- as.numeric(times)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  # Function to compute survival metric for a single time point
  predict_at_time <- function(t_val, th) {
    if (distribution == "exponential") {
      lambd <- th[1]
      gamma <- 1
    } else {
      lambd <- th[1]
      gamma <- th[2]
    }
    surv <- exp(-lambd * t_val^gamma)
    haz <- lambd * gamma * t_val^(gamma - 1)
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
  cov_m <- delta_method(theta, transform = g_func, covariance = covariance,
                        deriv_method = deriv_method, dx = dx)
  var_m <- diag(cov_m)

  # Confidence intervals. pmax(var_m, 0) guards the finite-difference paths,
  # where step-size cancellation can drive a near-zero variance slightly
  # negative; it is a no-op under exact autodiff, whose delta-method variance
  # is nonnegative by construction. Python does not clamp.
  se_m <- sqrt(pmax(var_m, 0))
  z_alpha <- qnorm(1 - alpha / 2)

  data.frame(
    time = times,
    predicted = est,
    variance = var_m,
    lower = est - z_alpha * se_m,
    upper = est + z_alpha * se_m
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
#' @export
aft_predictions_individual <- function(X, times, theta, distribution,
                                        measure = "survival") {
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

  # Linear predictor
  xbeta <- as.numeric(X %*% beta)

  # Compute survival metric at each time point
  results <- matrix(NA_real_, nrow = n, ncol = length(times))

  for (j in seq_along(times)) {
    t_val <- times[j]
    log_t <- log(t_val)
    eps <- (log_t - xbeta) / sigma
    hazard_scaler <- 1 / (sigma * t_val)

    if (distribution %in% c("exponential", "weibull")) {
      surv <- exp(-exp(eps))
      haz <- hazard_scaler * exp(eps)
    } else if (distribution %in% c("lognormal", "log-normal")) {
      surv <- 1 - pnorm(eps)
      haz <- hazard_scaler * dnorm(eps) / surv
    } else if (distribution %in% c("loglogistic", "log-logistic")) {
      surv <- 1 / (1 + exp(eps))
      haz <- hazard_scaler / (1 + exp(-eps))
    } else {
      cli::cli_abort("Distribution {.val {distribution}} is not supported.")
    }

    results[, j] <- convert_survival_measures(surv, haz, measure)
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
#'   `deriv_method = "exact"`. Default `1e-9`.
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
#' @examples
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
#' m <- MEstimator(stacked_equations = psi,
#'                 init = c(mean(log(t_obs)), 0, 0))
#' m <- estimate(m, solver = "nleqslv")
#'
#' aft_predictions_function(
#'   X = matrix(c(1, 1), nrow = 1), times = c(5, 10, 20),
#'   theta = m@theta, covariance = m@variance,
#'   distribution = "weibull", measure = "risk"
#' )
#'
#' @export
aft_predictions_function <- function(X, times, theta, covariance,
                                     distribution, measure = "survival",
                                     alpha = 0.05, deriv_method = "capprox",
                                     dx = 1e-9) {
  # Validate alpha level (mirrors the sibling prediction helpers)
  if (alpha <= 0 || alpha >= 1) {
    cli::cli_abort("{.arg alpha} must be between 0 and 1 (exclusive).")
  }

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

  # Measure for the single covariate pattern at one time point. Written with
  # scalar operations on the parameter vector so it differentiates exactly:
  # under exact autodiff `th` is a tangent-carrying pair vector, and indexing
  # plus scalar arithmetic keep the derivatives attached (whereas the numeric-
  # matrix path in aft_predictions_individual would strip them). The
  # autodiff-compatible standard_normal_cdf/standard_normal_pdf carry the
  # tangent through the log-normal error distribution. The same code returns
  # plain doubles when `th` is numeric, so the finite-difference paths reuse it.
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
      cli::cli_abort("Distribution {.val {distribution}} is not supported.")
    }
    convert_survival_measures(surv, haz, measure)
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
  cov_m <- delta_method(theta, transform = g_func, covariance = covariance,
                        deriv_method = deriv_method, dx = dx)
  var_m <- diag(as.matrix(cov_m))

  # Wald-type point-wise confidence intervals. pmax(var_m, 0) guards the finite-
  # difference paths, where step-size cancellation can drive a near-zero
  # variance slightly negative; it is a no-op under exact autodiff. Python does
  # not clamp.
  se_m <- sqrt(pmax(var_m, 0))
  z_alpha <- qnorm(1 - alpha / 2)

  data.frame(
    time = times,
    predicted = est,
    variance = var_m,
    lower = est - z_alpha * se_m,
    upper = est + z_alpha * se_m
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
#' @export
plogit_predict <- function(theta, time, event, X, S = NULL,
                            times_to_predict = NULL,
                            measure = "survival",
                            unique_times = NULL) {
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
  log_odds_w <- as.numeric(X %*% beta_x)            # n-vector
  log_odds_t <- as.numeric(time_design_matrix %*% beta_s)  # K-vector

  # Predicted event probability at each time: K-by-n matrix
  log_odds_w_matrix <- matrix(
    rep(log_odds_w, each = n_time_steps),
    nrow = n_time_steps, ncol = n
  )
  y_pred <- inverse_logit(log_odds_w_matrix + log_odds_t)

  # Cumulative product of (1 - p) gives survival: S(t) = prod(1 - h(t_k))
  survival_pred <- apply(1 - y_pred, 2, cumprod)     # K-by-n matrix

  # Convert to requested measure
  prediction_matrix <- convert_survival_measures(
    survival_pred, hazard = y_pred, measure = measure
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

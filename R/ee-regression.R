#' Estimating equation for regression
#'
#' Returns a p-by-n matrix of estimating equation contributions for
#' regression models. Supports linear, logistic, and Poisson regression:
#' \deqn{\psi_i(\theta) = \{Y_i - g(X_i^T \theta)\} X_i}
#'
#' @param theta Numeric vector of length p (number of covariates).
#' @param X Numeric n-by-p design matrix.
#' @param y Numeric vector of n observed outcome values.
#' @param model Character string: `"linear"`, `"logistic"`, or `"poisson"`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_regression,
#'   model = "linear"
#' )
#' summary(fit)
#'
#' # The same equation fits a logistic regression through the model argument.
#' fit_logit <- m_estimate(
#'   vs ~ mpg,
#'   data = mtcars,
#'   .ee = ee_regression,
#'   model = "logistic"
#' )
#' coef(fit_logit)
#'
#' @export
ee_regression <- function(theta, X, y, model, weights = NULL, offset = NULL) {
  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  w <- generate_weights(n, weights)

  # Linear predictor
  eta <- pt_as_vector(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }

  # Apply link function
  pred_y <- model_transform(eta, model)

  # Residuals and score
  residuals <- y - pred_y
  t(X * (w * residuals))
}

#' Estimating equation for ridge regression
#'
#' Returns a p-by-n matrix of estimating equation contributions for ridge
#' (L2-penalized) regression:
#' \deqn{\psi_i(\theta) = \{Y_i - g(X_i^T \theta)\} X_i - \frac{\lambda}{n} \theta}
#'
#' @param theta Numeric vector of length p.
#' @param X Numeric n-by-p design matrix.
#' @param y Numeric vector of n observed outcome values.
#' @param model Character string: `"linear"`, `"logistic"`, or `"poisson"`.
#' @param penalty Numeric scalar or vector of length p. Must be non-negative.
#'   Penalty terms scaled by n internally.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param center Numeric scalar or vector. Center for the penalty. Default `0`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' # A penalty vector gives one value per column of the design matrix. A scalar
#' # penalty would shrink the intercept along with the slopes.
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_ridge_regression,
#'   model = "linear",
#'   penalty = c(0, 5, 5)
#' )
#' coef(fit)
#'
#' @export
ee_ridge_regression <- function(
  theta,
  X,
  y,
  model,
  penalty,
  weights = NULL,
  center = 0,
  offset = NULL
) {
  # Ridge is the bridge penalty with gamma = 2 (L2 penalty)
  ee_bridge_regression(
    theta,
    X = X,
    y = y,
    model = model,
    weights = weights,
    penalty = penalty,
    gamma = 2,
    center = center,
    offset = offset
  )
}

#' Estimating equation for bridge penalized regression
#'
#' Returns a p-by-n matrix for bridge penalized regression. Bridge is the
#' general case: ridge is `gamma = 2`, LASSO approximation is
#' `gamma = 1 + epsilon`.
#'
#' @param theta Numeric vector of length p.
#' @param X Numeric n-by-p design matrix.
#' @param y Numeric vector of n observed outcome values.
#' @param model Character string: `"linear"`, `"logistic"`, or `"poisson"`.
#' @param penalty Numeric scalar or vector of length p. Must be non-negative.
#' @param gamma Numeric bridge exponent. Must be at least 1. Values below 2
#'   yield a penalty that is not everywhere differentiable, so the sandwich
#'   variance is not defined in all settings and a warning is issued.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param center Numeric scalar or vector. Default `0`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' # A penalty vector gives one value per column of the design matrix. A scalar
#' # penalty would shrink the intercept along with the slopes. The estimating
#' # equation carries the penalty's derivative, which varies like
#' # |theta|^(gamma - 1) near the penalty center. That derivative is itself
#' # differentiable only once gamma reaches 2, so gamma = 2.3 issues no
#' # warning while a value below 2 would.
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_bridge_regression,
#'   model = "linear",
#'   penalty = c(0, 5, 5),
#'   gamma = 2.3
#' )
#' coef(fit)
#'
#' @export
ee_bridge_regression <- function(
  theta,
  X,
  y,
  model,
  penalty,
  gamma,
  weights = NULL,
  center = 0,
  offset = NULL
) {
  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  eta <- pt_as_vector(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }
  pred_y <- model_transform(eta, model)

  # Unweighted p-by-n score, mirroring Python's ((y - pred_y) * X).T
  score <- t(X * (y - pred_y))

  # Length-p penalty recycles down the rows, subtracting from every column
  penalty_terms <- bridge_penalty(theta, n, penalty, gamma, center)
  psi <- score - penalty_terms

  # With no weights every observation contributes equally, so return the score
  # directly rather than multiplying by an n-by-p vector of ones. The result is
  # identical because multiplication by one is exact.
  if (is.null(weights)) {
    return(psi)
  }

  # Weight each observation's column so weights scale score and penalty
  # together. rep(w, each = ncol(X)) lays the length-n weights out column by
  # column to match the column-major p-by-n matrix, scaling every observation's
  # column in place rather than routing through two additional transposes.
  w <- generate_weights(n, weights)
  psi * rep(w, each = ncol(X))
}

#' Estimating equation for approximate LASSO regression
#'
#' Uses the bridge penalty with `gamma = 1 + epsilon` to approximate LASSO
#' (L1 penalty). The true LASSO is not differentiable at zero, so an
#' approximation is used.
#'
#' @inheritParams ee_bridge_regression
#' @param epsilon Numeric approximation parameter. Default `0.003`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' # A penalty vector gives one value per column of the design matrix. A scalar
#' # penalty would shrink the intercept along with the slopes.
#' #
#' # The approximate L1 penalty enters the estimating equation as its own
#' # derivative, and that derivative has unbounded slope at the penalty center.
#' # The estimating equation is therefore not differentiable there, so the
#' # bread matrix is undefined and every evaluation warns that the sandwich
#' # variance should not be trusted here.
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_lasso_regression,
#'   model = "linear",
#'   penalty = c(0, 5, 5)
#' )
#' coef(fit)
#'
#' @export
ee_lasso_regression <- function(
  theta,
  X,
  y,
  model,
  penalty,
  epsilon = 3e-3,
  weights = NULL,
  center = 0,
  offset = NULL
) {
  check_epsilon(epsilon)
  ee_bridge_regression(
    theta,
    X = X,
    y = y,
    model = model,
    weights = weights,
    penalty = penalty,
    gamma = 1 + epsilon,
    center = center,
    offset = offset
  )
}

#' Estimating equation for differentiable LASSO regression
#'
#' Uses a smooth approximation to the L1 penalty based on the standard
#' normal CDF and PDF.
#'
#' @inheritParams ee_bridge_regression
#' @param s Numeric smoothing parameter. Must be greater than zero.
#'   Default `1e-6`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' # A penalty vector gives one value per column of the design matrix. A scalar
#' # penalty would shrink the intercept along with the slopes. The estimating
#' # equation carries the penalty's derivative. Here that derivative is smooth,
#' # so the estimating equation is differentiable and no warning is issued,
#' # unlike ee_lasso_regression() whose derivative has unbounded slope at the
#' # penalty center.
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_dlasso_regression,
#'   model = "linear",
#'   penalty = c(0, 5, 5)
#' )
#' coef(fit)
#'
#' @export
ee_dlasso_regression <- function(
  theta,
  X,
  y,
  model,
  penalty,
  s = 1e-6,
  weights = NULL,
  center = 0,
  offset = NULL
) {
  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  eta <- pt_as_vector(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }
  pred_y <- model_transform(eta, model)

  # Unweighted p-by-n score, mirroring Python's ((y - pred_y) * X).T
  score <- t(X * (y - pred_y))

  # Length-p penalty recycles down the rows, subtracting from every column
  penalty_terms <- dlasso_penalty(theta, n, penalty, s, center)
  psi <- score - penalty_terms

  # With no weights every observation contributes equally, so return the score
  # directly rather than multiplying by an n-by-p vector of ones. The result is
  # identical because multiplication by one is exact.
  if (is.null(weights)) {
    return(psi)
  }

  # Weight each observation's column so weights scale score and penalty
  # together. rep(w, each = ncol(X)) lays the length-n weights out column by
  # column to match the column-major p-by-n matrix, scaling every observation's
  # column in place rather than routing through two additional transposes.
  w <- generate_weights(n, weights)
  psi * rep(w, each = ncol(X))
}

#' Estimating equation for elastic net regression
#'
#' Combines L1 (approximate LASSO via bridge) and L2 (ridge) penalties at
#' a given ratio. When `ratio = 1`, this is LASSO; when `ratio = 0`, ridge.
#'
#' @inheritParams ee_bridge_regression
#' @param ratio Numeric between 0 and 1. Proportion of L1 vs L2 penalty.
#' @param epsilon Numeric LASSO approximation parameter. Default `0.003`.
#'
#' @returns A p-by-n matrix.
#'
#' @examples
#' # A penalty vector gives one value per column of the design matrix. A scalar
#' # penalty would shrink the intercept along with the slopes.
#' #
#' # The L1 half of the penalty enters the estimating equation as its own
#' # derivative, and that derivative has unbounded slope at the penalty center.
#' # The estimating equation is therefore not differentiable there, so the
#' # bread matrix is undefined and every evaluation warns that the sandwich
#' # variance should not be trusted here.
#' fit <- m_estimate(
#'   mpg ~ wt + hp,
#'   data = mtcars,
#'   .ee = ee_elasticnet_regression,
#'   model = "linear",
#'   penalty = c(0, 5, 5),
#'   ratio = 0.5
#' )
#' coef(fit)
#'
#' @export
ee_elasticnet_regression <- function(
  theta,
  X,
  y,
  model,
  penalty,
  ratio,
  epsilon = 3e-3,
  weights = NULL,
  center = 0,
  offset = NULL
) {
  if (ratio < 0 || ratio > 1) {
    cli::cli_abort("The elastic-net {.arg ratio} must be between 0 and 1.")
  }
  check_epsilon(epsilon)

  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  eta <- pt_as_vector(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }
  pred_y <- model_transform(eta, model)

  # Unweighted p-by-n score, mirroring Python's ((y - pred_y) * X).T
  score <- t(X * (y - pred_y))

  # Combined penalty; length-p vector recycles down the rows
  penalty_l1 <- bridge_penalty(theta, n, penalty, gamma = 1 + epsilon, center)
  penalty_l2 <- bridge_penalty(theta, n, penalty, gamma = 2, center)
  penalty_terms <- ratio * penalty_l1 + (1 - ratio) * penalty_l2
  psi <- score - penalty_terms

  # With no weights every observation contributes equally, so return the score
  # directly rather than multiplying by an n-by-p vector of ones. The result is
  # identical because multiplication by one is exact.
  if (is.null(weights)) {
    return(psi)
  }

  # Weight each observation's column so weights scale score and penalty
  # together. rep(w, each = ncol(X)) lays the length-n weights out column by
  # column to match the column-major p-by-n matrix, scaling every observation's
  # column in place rather than routing through two additional transposes.
  w <- generate_weights(n, weights)
  psi * rep(w, each = ncol(X))
}

#' Internal bridge penalty calculation
#' @noRd
bridge_penalty <- function(theta, n, penalty, gamma, center) {
  penalty <- as.numeric(penalty)
  center <- as.numeric(center)

  check_penalty_shape(theta, penalty, center)

  if (gamma < 1) {
    cli::cli_abort(
      "{.code L_gamma} for {.arg gamma} < 1 cannot be supported with
       estimating equations evaluated using numerical methods."
    )
  }
  if (gamma < 2) {
    cli::cli_warn(
      "The estimating equation for the chosen penalized regression model is
       not always differentiable. Therefore, the bread matrix is not always
       defined for finite samples, and the sandwich should not be used to
       estimate the variance."
    )
  }

  penalty_scaled <- penalty / (gamma * n)
  penalty_scaled *
    gamma *
    (abs(theta - center)^(gamma - 1)) *
    sign(theta - center)
}

#' Internal differentiable LASSO penalty
#' @noRd
dlasso_penalty <- function(theta, n, penalty, s, center) {
  penalty <- as.numeric(penalty)
  center <- as.numeric(center)

  check_penalty_shape(theta, penalty, center)

  if (s <= 0) {
    cli::cli_abort(
      "{.arg s} must be greater than zero for the approximate LASSO."
    )
  }

  penalty_scaled <- penalty / n
  tc <- theta - center
  # The standard-normal wrappers carry tangents under exact differentiation and
  # reduce to pnorm/dnorm for plain numeric input, mirroring Python's use of
  # standard_normal_cdf and standard_normal_pdf here.
  penalty_scaled *
    (2 *
      standard_normal_cdf(tc / s) +
      2 * (tc / s) * standard_normal_pdf(tc / s) -
      1)
}

#' Estimating equation for multinomial logistic regression
#'
#' Supports unranked categorical outcomes. `y` must be an n-by-k indicator
#' matrix where the first column is the reference category. Returns a
#' `(b * (k-1))`-by-n matrix.
#'
#' @param theta Numeric vector of length `b * (k-1)`.
#' @param X Numeric n-by-b design matrix.
#' @param y Numeric n-by-k indicator matrix (first column = reference).
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A `(b * (k-1))`-by-n matrix.
#'
#' @examples
#' set.seed(123)
#' n <- 50
#' W <- rbinom(n, 1, 0.5)
#' probs <- cbind(0.5 - 0.2 * W, 0.3 + 0.1 * W, 0.2 + 0.1 * W)
#' y_cat <- sapply(seq_len(n), function(i) sample(1:3, 1, prob = probs[i, ]))
#'
#' # The outcome is an indicator matrix whose first column is the reference
#' # category.
#' y <- cbind(
#'   as.integer(y_cat == 1),
#'   as.integer(y_cat == 2),
#'   as.integer(y_cat == 3)
#' )
#' X <- cbind(1, W)
#'
#' psi <- function(theta) ee_mlogit(theta, X = X, y = y)
#'
#' # Two columns of X and two non-reference categories give four parameters.
#' m <- MEstimator(stacked_equations = psi, init = rep(0, 4)) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_mlogit <- function(theta, X, y, weights = NULL, offset = NULL) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  b <- ncol(X)
  k <- ncol(y)
  if (nrow(y) != n) {
    cli::cli_abort(
      "{.arg y} must have the same number of rows as {.arg X} ({n}), not
       {nrow(y)}."
    )
  }
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  w <- generate_weights(n, weights)

  off <- if (is.null(offset)) rep(0, n) else as.numeric(offset)

  expected_params <- b * (k - 1)
  if (length(theta) != expected_params) {
    cli::cli_abort(
      c(
        "Parameter length mismatch.",
        "x" = "Got {length(theta)} parameters, expected {expected_params}.",
        "i" = "For {k} categories and {b} predictors, need {b} * ({k} - 1) = {expected_params} parameters."
      )
    )
  }

  # Compute denominator: 1 + sum(exp(X %*% beta_j))
  denom <- rep(1, n)
  exp_pred <- vector("list", k - 1)
  for (j in seq_len(k - 1)) {
    idx <- ((j - 1) * b + 1):(j * b)
    beta_j <- theta[idx]
    exp_pred[[j]] <- exp(pt_as_vector(X %*% beta_j) + off)
    denom <- denom + exp_pred[[j]]
  }

  # Estimating equations for each non-reference category
  efuncs <- vector("list", k - 1)
  for (j in seq_len(k - 1)) {
    yhat_ref <- y[, 1] - 1 / denom
    y_j <- y[, j + 1]
    yhat_j <- yhat_ref + (y_j - exp_pred[[j]] / denom)
    efuncs[[j]] <- t(X * (w * yhat_j))
  }

  do.call(rbind, efuncs)
}

#' Estimating equation for beta regression
#'
#' Beta regression for outcomes in (0, 1) using mean-precision
#' parameterization. The last element of theta is log(phi), the log
#' precision parameter.
#'
#' @param theta Numeric vector of length `b + 1`.
#' @param X Numeric n-by-b design matrix.
#' @param y Numeric vector of n outcomes in (0, 1).
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A `(b+1)`-by-n matrix.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' set.seed(42)
#' n <- 50
#' W <- rnorm(n)
#' mu <- 1 / (1 + exp(-(0.5 + 0.3 * W)))
#' y <- rbeta(n, shape1 = mu * 10, shape2 = (1 - mu) * 10)
#' X <- cbind(1, W)
#'
#' psi <- function(theta) ee_beta_regression(theta, X = X, y = y)
#'
#' # The last parameter is log(phi), started here at a precision of 10. The
#' # default solver diverges on this equation, so use nleqslv.
#' m <- MEstimator(stacked_equations = psi, init = c(0, 0, log(10))) |>
#'   estimate(solver = "nleqslv")
#' coef(m)
#'
#' @export
ee_beta_regression <- function(theta, X, y, weights = NULL, offset = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  b <- ncol(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  w <- generate_weights(n, weights)

  off <- if (is.null(offset)) rep(0, n) else as.numeric(offset)

  beta <- theta[seq_len(b)]
  phi <- exp(theta[b + 1])

  yhat <- inverse_logit(pt_as_vector(X %*% beta) + off)
  # Clamp yhat away from exact 0/1 to avoid digamma poles during iteration
  yhat <- pmin(pmax(yhat, 1e-15), 1 - 1e-15)
  logit_y <- logit(y)

  resid <- logit_y - deli_digamma(yhat * phi) + deli_digamma((1 - yhat) * phi)

  # Regression coefficients EE
  ef_mean <- t(X * (w * yhat * (1 - yhat) * phi * resid))

  # Precision parameter EE
  ef_prc <- matrix(
    w *
      (deli_digamma(phi) -
        yhat * deli_digamma(yhat * phi) -
        (1 - yhat) * deli_digamma((1 - yhat) * phi) +
        yhat * log(y) +
        (1 - yhat) * log(1 - y)),
    nrow = 1
  )

  rbind(ef_mean, ef_prc)
}

#' Estimating equation for Tobit regression (Type I)
#'
#' Handles left and/or right censored outcomes using standard normal PDF/CDF.
#' Theta is `(beta, log(sigma))`, a vector of length `b + 1`.
#'
#' @param theta Numeric vector of length `b + 1`.
#' @param X Numeric n-by-b design matrix.
#' @param y Numeric vector of n observed (possibly censored) outcome values.
#' @param lower Numeric lower censoring limit, or `NULL` (default, no left
#'   censoring).
#' @param upper Numeric upper censoring limit, or `NULL` (default, no right
#'   censoring).
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A `(b+1)`-by-n matrix.
#'
#' @examples
#' # A latent outcome observed only down to zero, so the negative values are
#' # left censored at the limit.
#' set.seed(123)
#' n <- 200
#' X <- cbind(1, rnorm(n))
#' y <- pmax(1 + 0.5 * X[, 2] + rnorm(n), 0)
#'
#' psi <- function(theta) ee_tobit(theta, X = X, y = y, lower = 0)
#'
#' # The last parameter is log(sigma), started at the log of the observed
#' # standard deviation.
#' m <- MEstimator(
#'   stacked_equations = psi,
#'   init = c(mean(y), 0, log(sd(y)))
#' ) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_tobit <- function(
  theta,
  X,
  y,
  lower = NULL,
  upper = NULL,
  weights = NULL,
  offset = NULL
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  b <- ncol(X)
  check_data_length(y, n, "y")
  if (!is.null(offset)) {
    check_data_length(offset, n, "offset")
  }
  w <- generate_weights(n, weights)

  beta <- theta[seq_len(b)]
  sigma <- exp(theta[b + 1])

  # Linear predictor
  yhat <- pt_as_vector(X %*% beta)
  if (!is.null(offset)) {
    yhat <- yhat + as.numeric(offset)
  }
  resid <- y - yhat

  # Left censoring
  if (!is.null(lower)) {
    lcensor <- as.numeric(y <= lower)
    scaled_yl <- (lower - yhat) / sigma
    pdf_lower <- standard_normal_pdf(scaled_yl)
    cdf_lower <- pmax(standard_normal_cdf(scaled_yl), 1e-14)
    lambda_lower <- pdf_lower / cdf_lower
  } else {
    lower <- -Inf
    lcensor <- 0
    lambda_lower <- 0
    scaled_yl <- 1
  }

  # Right censoring
  if (!is.null(upper)) {
    rcensor <- as.numeric(y >= upper)
    scaled_yu <- (upper - yhat) / sigma
    pdf_upper <- standard_normal_pdf(scaled_yu)
    cdf_upper <- pmax(1 - standard_normal_cdf(scaled_yu), 1e-14)
    lambda_upper <- pdf_upper / cdf_upper
  } else {
    upper <- Inf
    rcensor <- 0
    lambda_upper <- 0
    scaled_yu <- 1
  }

  ucensor <- (1 - lcensor) * (1 - rcensor)

  # Input validation
  if (lower >= upper) {
    cli::cli_abort(
      "The {.arg lower} limit must be less than the {.arg upper} limit."
    )
  }
  if (any(y < lower)) {
    cli::cli_abort(
      "Some observations have values below the specified {.arg lower} limit."
    )
  }
  if (any(y > upper)) {
    cli::cli_abort(
      "Some observations have values above the specified {.arg upper} limit."
    )
  }

  # Regression score: (b)-by-n
  ef_treg <- t(
    X *
      (w *
        (lcensor *
          (-lambda_lower / sigma) +
          ucensor * (resid / sigma^2) +
          rcensor * (lambda_upper / sigma)))
  )

  # Variance score: 1-by-n
  ef_sigma <- matrix(
    w *
      (lcensor *
        (-scaled_yl * lambda_lower / sigma) +
        ucensor * (-1 / sigma + resid^2 / sigma^3) +
        rcensor * (scaled_yu * lambda_upper / sigma)),
    nrow = 1
  )

  rbind(ef_treg, ef_sigma)
}

#' Estimating equation for additive regression (GAM)
#'
#' Generalized Additive Model via L2-penalized splines. Internally expands `X`
#' using [additive_design_matrix()] and delegates to [ee_bridge_regression()]
#' with `gamma = 2` (ridge penalty). The penalty only applies to the spline
#' basis terms, not to the original linear terms.
#'
#' @param theta Numeric vector of length equal to the number of columns in
#'   the expanded additive design matrix.
#' @param X Numeric n-by-b design matrix (before spline expansion).
#' @param y Numeric vector of n observed outcome values.
#' @param specifications A list of length `b` controlling spline generation.
#'   Each element is either `NULL` (no spline) or a list with keys `knots`,
#'   and optionally `natural`, `power`, `penalty`, `normalized`. See
#'   [additive_design_matrix()] for details.
#' @param model Character string: `"linear"`, `"logistic"`, or `"poisson"`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix, where p is the number of columns in the expanded
#'   additive design matrix.
#'
#' @examples
#' set.seed(42)
#' n <- 200
#' x <- runif(n, -3, 3)
#' y <- sin(x) + rnorm(n, sd = 0.3)
#' X <- cbind(1, x)
#'
#' # No spline on the intercept column, a penalized spline on x.
#' specs <- list(NULL, list(knots = c(-2, -1, 0, 1, 2), penalty = 5))
#'
#' psi <- function(theta) {
#'   ee_additive_regression(
#'     theta,
#'     X = X,
#'     y = y,
#'     specifications = specs,
#'     model = "linear"
#'   )
#' }
#'
#' # One parameter per column of the expanded design matrix.
#' m <- MEstimator(
#'   stacked_equations = psi,
#'   init = rep(0, ncol(additive_design_matrix(X, specs)))
#' ) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_additive_regression <- function(
  theta,
  X,
  y,
  specifications,
  model,
  weights = NULL,
  offset = NULL
) {
  # Build the additive design matrix and retrieve per-column penalty values
  result <- additive_design_matrix(
    X = X, # Original design matrix
    specifications = specifications, # Spline specs per column
    return_penalty = TRUE
  ) # Also return penalty vector

  Xa <- result$X # Expanded design matrix with spline columns
  penalty <- result$penalty # Penalty vector (0 for linear terms, lambda for spline terms)

  # Delegate to bridge regression with gamma = 2 (L2 / ridge penalty)
  ee_bridge_regression(
    theta = theta,
    X = Xa, # Use the expanded design matrix
    y = y, # Observed outcomes
    model = model, # Link function
    penalty = penalty, # Per-column penalty values
    gamma = 2, # Ridge (L2) penalty
    weights = weights, # Observation weights
    center = 0, # Splines always penalized toward zero
    offset = offset
  ) # Optional offset
}

#' Internal model transform dispatcher
#' @noRd
model_transform <- function(eta, model) {
  model <- tolower(model)
  if (model == "linear") {
    eta
  } else if (model == "logistic") {
    1 / (1 + exp(-eta))
  } else if (model == "poisson") {
    exp(eta)
  } else {
    cli::cli_abort(
      "Model {.val {model}} is not supported. Use {.val linear}, {.val logistic}, or {.val poisson}."
    )
  }
}

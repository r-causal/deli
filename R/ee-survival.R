#' Estimating equation for accelerated failure time models
#'
#' Returns a p-by-n matrix of estimating equation contributions for
#' accelerated failure time (AFT) models. Supports Weibull, exponential,
#' log-logistic, and log-normal distributions.
#'
#' The AFT model uses the parameterization where
#' \eqn{Z_i = (\log(t_i) - X_i \beta) / \sigma} and the last element of
#' theta is \eqn{\log(1/\sigma)}.
#'
#' @param theta Numeric vector of length `b + 1` (or `b` for exponential).
#'   The first `b` elements are regression coefficients and the last element
#'   is `log(1/sigma)` (the log-inverse scale).
#' @param X Numeric n-by-b design matrix (should include intercept column).
#' @param time Numeric vector of n observed (possibly censored) times.
#' @param event Numeric vector of n event indicators (1 = event, 0 = censored).
#' @param distribution Character string: `"weibull"`, `"exponential"`,
#'   `"log-logistic"`, or `"log-normal"`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix where p is the number of parameters.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' # A Weibull AFT model for times generated from one binary covariate, with
#' # some of the times right censored by an independent exponential.
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
#'
#' # The default rootSolve solver does not converge here, so nleqslv is used.
#' # The last parameter is log(1/sigma), not a regression coefficient.
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(mean(log(t_obs)), 0, 0),
#'   solver = "nleqslv"
#' )
#' coef(m)
#'
#' @export
ee_aft <- function(
  theta,
  X,
  time,
  event,
  distribution,
  weights = NULL,
  offset = NULL
) {
  X <- as.matrix(X)
  time <- as.numeric(time)
  event <- as.numeric(event)
  n <- nrow(X)
  beta_dim <- ncol(X)
  distribution <- tolower(distribution)

  # Validate survival data
  check_survival_data_valid(delta = event, time = time)

  # Extract coefficients
  beta <- theta[seq_len(beta_dim)]

  # Sigma: for exponential, sigma is fixed at 1

  if (distribution == "exponential") {
    sigma <- 1
  } else {
    sigma <- exp(-theta[length(theta)])
  }

  # Linear predictor with optional offset
  xbeta <- pt_as_vector(X %*% beta)
  if (!is.null(offset)) {
    xbeta <- xbeta + as.numeric(offset)
  }

  # Computing error distribution for each observation
  z_i <- (log(time) - xbeta) / sigma

  # Weights
  w <- generate_weights(n, weights)

  # Handling different distribution specifications
  if (distribution %in% c("exponential", "weibull")) {
    df_f <- 1 - exp(z_i) # f'/f for Weibull/exponential
    dS_S <- exp(z_i) # S'/S for Weibull/exponential
  } else if (distribution %in% c("log-logistic", "loglogistic")) {
    df_f <- 1 - (2 * exp(z_i)) / (1 + exp(z_i)) # f'/f for log-logistic
    dS_S <- exp(z_i) / (1 + exp(z_i)) # S'/S for log-logistic
  } else if (distribution %in% c("log-normal", "lognormal")) {
    df_f <- -z_i # f'/f for log-normal
    dS_S <- standard_normal_pdf(z_i) / (1 - standard_normal_cdf(z_i)) # S'/S for log-normal
  } else {
    cli::cli_abort(
      c(
        "The distribution {.val {distribution}} is not supported.",
        "i" = "Use one of: {.val exponential}, {.val weibull},
               {.val log-logistic}, {.val log-normal}."
      )
    )
  }

  # Individual contributions
  lambda_epsilon <- event * df_f - (1 - event) * dS_S

  # Score for regression coefficients: (-1/sigma) * lambda_epsilon * X
  score_scale <- (-1 / sigma) * lambda_epsilon # n-length vector
  ef_beta <- t(X * (w * score_scale)) # b-by-n matrix

  # Return based on distribution

  if (distribution == "exponential") {
    return(ef_beta)
  }

  # Score for shape parameter: (-1/sigma) * lambda_epsilon * z_i - event/sigma
  score_shape <- (-1 / sigma) * lambda_epsilon * z_i - event / sigma
  ef_shape <- matrix(w * score_shape, nrow = 1) # 1-by-n matrix

  rbind(ef_beta, ef_shape)
}

#' Estimating equation for parametric survival models
#'
#' Returns a p-by-n matrix of estimating equation contributions for parametric
#' survival models. Supports exponential, Weibull, and Gompertz distributions.
#'
#' The estimating equations are based on the score equations of the
#' corresponding parametric survival model, accounting for right censoring.
#' For event observations, the contribution comes from the log-density; for
#' censored observations, the contribution comes from the log-survival function.
#'
#' @param theta Numeric vector of distribution parameters. For exponential,
#'   a single parameter (lambda). For Weibull and Gompertz, two parameters
#'   (lambda, gamma).
#' @param time Numeric vector of n observed (possibly censored) times.
#' @param event Numeric vector of n event indicators (1 = event, 0 = censored).
#' @param distribution Character string: `"exponential"`, `"weibull"`, or
#'   `"gompertz"`.
#'
#' @returns A p-by-n matrix where p is the number of parameters.
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
#' # The first parameter is the scale, the second the shape. A shape near 1
#' # means the Weibull fits about as well as the simpler exponential.
#' coef(m)
#'
#' @export
ee_survival_model <- function(theta, time, event, distribution) {
  time <- as.numeric(time)
  event <- as.numeric(event)
  distribution <- tolower(distribution)

  # Validate survival data
  check_survival_data_valid(delta = event, time = time)

  # Extract parameters based on distribution
  if (distribution == "exponential") {
    lambd <- theta[1]
    gamma <- 1
  } else {
    lambd <- theta[1]
    gamma <- theta[2]
  }

  # Compute estimating equation contributions by distribution
  if (distribution %in% c("exponential", "weibull")) {
    # Estimating equation for lambda
    ef_lambda <- (event / lambd) - time^gamma
    # Estimating equation for gamma
    ef_gamma <- (event / gamma) +
      (event * log(time)) -
      (lambd * (time^gamma) * log(time))
  } else if (distribution == "gompertz") {
    exp_gt <- exp(gamma * time) # exp(gamma * time)
    ef_lambda <- event / lambd - (exp_gt - 1) / gamma
    ef_gamma <- lambd /
      (gamma^2) *
      (exp_gt - 1) +
      event * time -
      (lambd / gamma) * exp_gt * time
  } else {
    cli::cli_abort(
      c(
        "The distribution {.val {distribution}} is not supported.",
        "i" = "Use one of: {.val exponential}, {.val weibull}, {.val gompertz}."
      )
    )
  }

  # Return stacked estimating equations
  if (distribution == "exponential") {
    return(matrix(ef_lambda, nrow = 1))
  }

  rbind(
    matrix(ef_lambda, nrow = 1),
    matrix(ef_gamma, nrow = 1)
  )
}

#' Estimating equation for pooled logistic regression
#'
#' Returns a p-by-n matrix of estimating equation contributions for pooled
#' logistic regression with discrete-time survival data. This implementation
#' does not require creation of a long data set.
#'
#' @param theta Numeric vector of length `b + p_s`, where `b` is the number
#'   of covariate columns in `X`. When `S` is supplied, `p_s` is `ncol(S)`;
#'   when `S = NULL`, `p_s` is `K`, the number of unique event times.
#' @param X Numeric n-by-b design matrix for baseline covariates.
#' @param time Numeric vector of n observed (possibly censored) times.
#' @param event Numeric vector of n event indicators (1 = event, 0 = censored).
#' @param S Optional time design matrix with K rows (one per time step) and
#'   p_s columns. Default `NULL` uses disjoint indicators for unique event
#'   times.
#' @param unique_times Optional numeric vector of unique event times. Only
#'   used when `S = NULL`. Default `NULL`.
#' @param weights Optional numeric vector of n weights, or an n-by-K matrix of
#'   time-varying weights with one column per time interval (K must equal the
#'   number of unit-time intervals). Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A `(b + p_s)`-by-n matrix.
#'
#' @examples
#' # Bladder tumor recurrence, comparing the novel treatment to placebo while
#' # adjusting for the number and size of the initial tumors.
#' d <- collett_bladder
#' d$novel <- d$treat - 1
#' W <- as.matrix(d[, c("novel", "init", "size")])
#'
#' # Time is modeled with disjoint indicators, one per distinct event time,
#' # which ee_plogit builds by default.
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
#' # The first three parameters are the covariate coefficients, which
#' # approximate log hazard ratios. The remaining parameters are the
#' # pooled-logistic time intercepts.
#' summary(m, subset = 1:3)
#'
#' @export
ee_plogit <- function(
  theta,
  X,
  time,
  event,
  S = NULL,
  unique_times = NULL,
  weights = NULL,
  offset = NULL
) {
  X <- as.matrix(X)
  time <- as.numeric(time) # observed times
  event <- as.numeric(event) # event indicator
  n <- nrow(X)
  xp <- ncol(X) # number of covariate parameters

  # Split theta into covariate and time parameters
  beta_x <- theta[seq_len(xp)]
  beta_s <- theta[(xp + 1):length(theta)]

  # Processing design matrix for time
  if (is.null(S)) {
    # Build disjoint indicator matrix from unique event times
    if (is.null(unique_times)) {
      event_times <- time[event == 1] # look up event times
      unique_times <- sort(unique(event_times)) # extract unique ones
    } else {
      unique_times <- as.numeric(unique_times)
    }
    n_time_steps <- length(unique_times)
    # Create disjoint indicator design matrix
    time_design_matrix <- diag(n_time_steps)
    time_design_matrix[, 1] <- 1 # first column is intercept
  } else {
    S <- as.matrix(S)
    time_design_matrix <- S
    unique_times <- seq_len(max(time)) # 1:max(time)
    n_time_steps <- length(unique_times)
    if (n_time_steps != nrow(time_design_matrix)) {
      cli::cli_abort(
        c(
          "Dimension mismatch between time intervals and {.arg S}.",
          "x" = "Found {n_time_steps} unit-time intervals but {.arg S} has
                 {nrow(time_design_matrix)} rows.",
          "i" = "These values must match."
        )
      )
    }
  }

  # Log-odds contributions for covariates and time
  log_odds_w <- pt_as_vector(X %*% beta_x) # n-length vector
  if (!is.null(offset)) {
    log_odds_w <- log_odds_w + as.numeric(offset)
  }
  log_odds_t <- pt_as_vector(time_design_matrix %*% beta_s) # K-length vector

  # Computing residuals across time intervals
  # log_odds_w_matrix: K-by-n matrix (stacked copies of covariate contributions)
  log_odds_w_matrix <- matrix(
    rep(log_odds_w, each = n_time_steps),
    nrow = n_time_steps,
    ncol = n
  )

  # y_obs: K-by-n indicator matrix (event at specific time interval)
  # Each row k: event_i * I(time_i == unique_times[k])
  # Use rep(event, each=) to broadcast the indicator across rows (as for weights)
  y_obs <- rep(event, each = n_time_steps) * outer(unique_times, time, "==")

  # y_pred: K-by-n predicted probability matrix
  y_pred <- inverse_logit(log_odds_w_matrix + log_odds_t)

  # in_risk_set: K-by-n indicator (individual is at risk at time k)
  in_risk_set <- outer(unique_times, time, "<=")
  storage.mode(in_risk_set) <- "double"

  # Residual matrix: K-by-n
  residual_matrix <- (y_obs - y_pred) * in_risk_set

  # Incorporating weights
  w <- generate_weights(n, weights)
  # Route every matrix through the column-count check, matching Python, which
  # rejects any 2D weight array whose column count differs from the number of
  # time intervals. An n-by-1 matrix is not treated as a vector: with more than
  # one interval it fails the check, as a genuine time-varying matrix would.
  if (is.matrix(w)) {
    # Time-varying weights: n-by-K -> transpose to K-by-n
    if (ncol(w) != n_time_steps || nrow(w) != n) {
      cli::cli_abort(
        c(
          "2D {.arg weights} dimension mismatch.",
          "x" = "Expected {n} row{?s} and {n_time_steps} column{?s},
                 got {nrow(w)} row{?s} and {ncol(w)} column{?s}."
        )
      )
    }
    w_matrix <- t(w) # K-by-n
    residual_matrix <- residual_matrix * w_matrix
  } else {
    # 1D weights: broadcast across time intervals
    residual_matrix <- residual_matrix * rep(w, each = n_time_steps)
  }

  # Score for X: sum residuals across time intervals, then multiply by X
  y_resid <- colSums(residual_matrix) # n-length vector
  x_score <- t(X * y_resid) # b-by-n matrix

  # Score for S (time parameters)
  if (is.null(S)) {
    # Disjoint indicators: residual_matrix is already K-by-n
    t_score <- residual_matrix
  } else {
    # Matrix multiplication: p_s-by-K %*% K-by-n -> p_s-by-n
    t_score <- t(S) %*% residual_matrix
  }

  # Stack covariate and time scores
  rbind(x_score, t_score)
}

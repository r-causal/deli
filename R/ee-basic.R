#' Estimating equation for the mean
#'
#' Returns a 1-by-n matrix of estimating equation contributions for the mean:
#' \eqn{\psi_i(\theta) = w_i (Y_i - \theta)}.
#'
#' @param theta Numeric vector of length 1.
#' @param y Numeric vector of observed values.
#' @param weights Optional numeric vector of weights (same length as `y`).
#'   Default `NULL` assigns weight 1 to all observations.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
#' psi <- function(theta) ee_mean(theta, y = y)
#' m <- MEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_mean <- function(theta, y, weights = NULL) {
  y <- as.numeric(y)
  w <- generate_weights(length(y), weights)
  matrix(w * (y - theta[1]), nrow = 1)
}

#' Estimating equations for the mean and variance
#'
#' Returns a 2-by-n matrix of estimating equation contributions for the mean
#' and variance:
#' \deqn{\psi_i(\theta) = \begin{pmatrix} Y_i - \theta_1 \\ (Y_i - \theta_1)^2 - \theta_2 \end{pmatrix}}
#'
#' @param theta Numeric vector of length 2. `theta[1]` is the mean,
#'   `theta[2]` is the variance.
#' @param y Numeric vector of observed values.
#'
#' @returns A 2-by-n matrix.
#'
#' @examples
#' y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
#' psi <- function(theta) ee_mean_variance(theta, y = y)
#' m <- MEstimator(stacked_equations = psi, init = c(0, 1)) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_mean_variance <- function(theta, y) {
  y <- as.numeric(y)
  rbind(
    y - theta[1],
    (y - theta[1])^2 - theta[2]
  )
}

#' Estimating equation for the robust mean
#'
#' Returns a 1-by-n matrix of estimating equation contributions for the robust
#' mean using the specified loss function:
#' \eqn{\psi_i(\theta) = f_k(Y_i - \theta)}.
#'
#' @param theta Numeric vector of length 1.
#' @param y Numeric vector of observed values.
#' @param k Numeric tuning parameter for the loss function.
#' @param loss Character string specifying the loss function. Default
#'   `"huber"`. See [robust_loss_functions()] for options.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' # Forty-nine standard normal observations plus one gross outlier.
#' set.seed(1)
#' y <- c(rnorm(49), 50)
#' psi <- function(theta) ee_mean_robust(theta, y = y, k = 1.345, loss = "huber")
#' m <- MEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#'
#' # The Huber estimate stays near the uncontaminated mean, unlike mean(y).
#' coef(m)
#' mean(y)
#'
#' @export
ee_mean_robust <- function(theta, y, k, loss = "huber") {
  y <- as.numeric(y)
  matrix(robust_loss_functions(y - theta[1], loss = loss, k = k), nrow = 1)
}

#' Estimating equation for the geometric mean
#'
#' Returns a 1-by-n matrix of estimating equation contributions for the
#' geometric mean. When `log_theta = TRUE` (default), solves
#' \eqn{\log(Y_i) - \log(\theta)}. When `log_theta = FALSE`, solves
#' \eqn{\log(Y_i) - \theta} where \eqn{\theta} is the log of the geometric
#' mean.
#'
#' @param theta Numeric vector of length 1.
#' @param y Numeric vector of positive observed values.
#' @param weights Optional numeric vector of weights. Default `NULL`.
#' @param log_theta Logical. If `TRUE` (default), internally log-transforms
#'   theta. Default `TRUE`.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
#' psi <- function(theta) ee_mean_geometric(theta, y = y)
#' m <- MEstimator(stacked_equations = psi, init = 1) |>
#'   estimate()
#' coef(m)
#'
#' @export
ee_mean_geometric <- function(theta, y, weights = NULL, log_theta = TRUE) {
  y <- as.numeric(y)
  w <- generate_weights(length(y), weights)
  if (log_theta) {
    matrix(w * (log(y) - log(theta[1])), nrow = 1)
  } else {
    matrix(w * (log(y) - theta[1]), nrow = 1)
  }
}

#' Estimating equation for the percentile
#'
#' Returns a 1-by-n matrix for the q-th percentile:
#' \eqn{\psi_i(\theta) = q - I(Y_i \le \theta)}.
#'
#' @details The derivative of this estimating equation is not defined at
#'   \eqn{\hat{\theta}}, so the bread matrix and sandwich variance cannot be
#'   used to estimate the variance. The function warns on each call for this
#'   reason. It is offered for completeness but is not generally recommended
#'   for applications.
#'
#' @param theta Numeric vector of length 1.
#' @param y Numeric vector of observed values.
#' @param q Numeric percentile, must be in `(0, 1)`.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
#' psi <- function(theta) ee_percentile(theta, y = y, q = 0.5)
#'
#' # The estimating function is a step function of theta, so its derivative is
#' # zero wherever it is defined and the root finder cannot search. Start at the
#' # sample quantile, which is the solution.
#' m <- MEstimator(stacked_equations = psi, init = median(y)) |>
#'   estimate()
#'
#' # Each call warns for the same reason: the median estimating equation is not
#' # differentiable, so the sandwich variance should not be trusted.
#' coef(m)
#'
#' @export
ee_percentile <- function(theta, y, q) {
  y <- as.numeric(y)
  if (q <= 0 || q >= 1) {
    cli::cli_abort("{.arg q} must be between 0 and 1 (exclusive).")
  }
  cli::cli_warn(
    "The estimating equation is not differentiable at {.arg theta}. Therefore,
     the bread matrix is not defined for finite samples, and the sandwich should
     not be used to estimate the variance."
  )
  matrix(q - as.numeric(y <= theta[1]), nrow = 1)
}

#' Estimating equations for the positive mean deviation
#'
#' Returns a 2-by-n matrix for the positive mean deviation and median.
#'
#' @details The derivative of the estimating equation for the median is not
#'   defined at \eqn{\hat{\theta}}, so the bread matrix and sandwich variance
#'   cannot be used to estimate the variance. The function warns on each call
#'   for this reason. It is offered for completeness but is not generally
#'   recommended for applications.
#'
#' @param theta Numeric vector of length 2. `theta[1]` is the positive mean
#'   deviation, `theta[2]` is the median.
#' @param y Numeric vector of observed values.
#'
#' @returns A 2-by-n matrix.
#'
#' @examplesIf requireNamespace("minpack.lm", quietly = TRUE)
#' y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
#' psi <- function(theta) ee_positive_mean_deviation(theta, y = y)
#'
#' # The median equation is a step function of theta, so the Jacobian is
#' # singular and the default solver cannot move. Start the median at the sample
#' # median, which is the solution, and use the Levenberg-Marquardt solver to
#' # fit the remaining parameter.
#' m <- MEstimator(stacked_equations = psi, init = c(0, median(y))) |>
#'   estimate(solver = "lm")
#'
#' # Each call warns for the same reason: the median estimating equation is not
#' # differentiable, so the sandwich variance should not be trusted.
#' coef(m)
#'
#' @export
ee_positive_mean_deviation <- function(theta, y) {
  y <- as.numeric(y)
  cli::cli_warn(
    "The estimating equation for the median is not differentiable. Therefore,
     the bread matrix is not defined for finite samples, and the sandwich should
     not be used to estimate the variance."
  )
  rbind(
    2 * (y - theta[2]) * as.numeric(y > theta[2]) - theta[1],
    0.5 - as.numeric(y <= theta[2])
  )
}

#' Estimating equation for E-max dose-response model
#'
#' Implements the hyperbolic E-max (Hill) model:
#' \deqn{R_i = \theta_0 + \frac{\theta_{max} D_i}{\theta_{50} + D_i}}
#'
#' @param theta Numeric vector of length 3: zero-dose response (`e0`),
#'   maximum response (`emax`), ED50.
#' @param dose Numeric vector of n dose values.
#' @param response Numeric vector of n response values.
#' @param loss Optional character string for robust loss function. Default
#'   `NULL` (no robust loss). See [robust_loss_functions()].
#' @param k Optional numeric tuning parameter for robust loss.
#'
#' @returns A 3-by-n matrix.
#'
#' @examples
#' # Dose-response of a herbicide on ryegrass root length. The response falls
#' # with dose, so the maximum response parameter is initialized negative.
#' psi <- function(theta) {
#'   ee_emax(theta, dose = inderjit$dose, response = inderjit$response)
#' }
#'
#' m <- m_estimate(stacked_equations = psi, init = c(8, -8, 2))
#'
#' # Zero-dose response, maximum change in response, and ED50
#' coef(m)
#'
#' @export
ee_emax <- function(theta, dose, response, loss = NULL, k = NULL) {
  dose <- as.numeric(dose)
  response <- as.numeric(response)

  e0 <- theta[1]
  emax <- theta[2]
  e50 <- theta[3]

  # Predicted response
  yhat <- e0 + (emax * dose) / (e50 + dose)
  resid <- response - yhat

  # Apply robust loss if specified
  if (!is.null(loss)) {
    if (is.null(k)) {
      cli::cli_abort(
        "{.arg k} must be provided when using a robust {.arg loss}."
      )
    }
    resid <- robust_loss_functions(resid, loss = loss, k = k)
  }

  # Gradient components
  d_e0 <- rep(1, length(dose))
  d_emax <- dose / (e50 + dose)
  d_e50 <- -(emax * dose) / (e50 + dose)^2

  rbind(
    matrix(resid * d_e0, nrow = 1),
    matrix(resid * d_emax, nrow = 1),
    matrix(resid * d_e50, nrow = 1)
  )
}

#' Estimating equation for delta-effective dose (E-max)
#'
#' Computes the effective dose at level delta from the E-max model.
#' Should be stacked with [ee_emax()].
#'
#' @param theta Numeric scalar. The ED(delta) parameter.
#' @param dose Numeric vector (used for dimension only).
#' @param delta Numeric effective dose level of interest.
#' @param ed50 Numeric estimated ED50 from E-max model.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' # This equation carries no information on its own, so stack it with the
#' # E-max equation that supplies the ED50. Here delta = 0.9 requests the ED90.
#' psi <- function(theta) {
#'   emax <- ee_emax(
#'     theta[1:3],
#'     dose = inderjit$dose,
#'     response = inderjit$response
#'   )
#'   ed90 <- ee_emax_ed(
#'     theta[4],
#'     dose = inderjit$dose,
#'     delta = 0.9,
#'     ed50 = theta[3]
#'   )
#'   rbind(emax, ed90)
#' }
#'
#' m <- m_estimate(stacked_equations = psi, init = c(8, -8, 2, 10))
#'
#' # theta_4 is the ED90, and stacking gives it its own sandwich standard error
#' coef(m)
#'
#' @export
ee_emax_ed <- function(theta, dose, delta, ed50) {
  n <- length(dose)
  ed_delta <- delta / (1 - delta) * ed50 - theta[1]
  matrix(rep(ed_delta, n), nrow = 1)
}

#' Estimating equation for 4-parameter log-logistic dose-response model
#'
#' Implements the 4-parameter log-logistic model:
#' \deqn{R_i = \theta_0 + \frac{\theta_m - \theta_0}{1 + \exp[\theta_s (\log D_i - \log \theta_{50})]}}
#'
#' @param theta Numeric vector of length 4: lower limit, upper limit,
#'   ED50, steepness.
#' @param dose Numeric vector of n dose values.
#' @param response Numeric vector of n response values.
#' @param loss Optional character string for robust loss function.
#' @param k Optional numeric tuning parameter for robust loss.
#'
#' @returns A 4-by-n matrix.
#'
#' @examplesIf requireNamespace("nleqslv", quietly = TRUE)
#' # Four-parameter log-logistic dose-response for the ryegrass data.
#' psi <- function(theta) {
#'   ee_loglogistic(theta, dose = inderjit$dose, response = inderjit$response)
#' }
#'
#' # The default rootSolve solver does not converge here, so nleqslv is used.
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(0.2, 8, 2, 1),
#'   solver = "nleqslv"
#' )
#'
#' # Lower limit, upper limit, ED50, and steepness
#' coef(m)
#'
#' @export
ee_loglogistic <- function(theta, dose, response, loss = NULL, k = NULL) {
  dose <- as.numeric(dose)
  response <- as.numeric(response)

  e0 <- theta[1] # lower limit
  emax <- theta[2] # upper limit
  e50 <- theta[3] # ED50
  steep <- theta[4] # steepness

  # Replace zero doses with small value to avoid log(0)
  dose <- ifelse(dose == 0, 1e-9, dose)

  rho <- (dose / e50)^steep
  yhat <- e0 + (emax - e0) / (1 + rho)
  resid <- response - yhat

  # Apply robust loss if specified
  if (!is.null(loss)) {
    if (is.null(k)) {
      cli::cli_abort(
        "{.arg k} must be provided when using a robust {.arg loss}."
      )
    }
    resid <- robust_loss_functions(resid, loss = loss, k = k)
  }

  # Gradient components
  d_lower <- 1 - 1 / (1 + rho)
  d_upper <- 1 / (1 + rho)
  d_ed50 <- (emax - e0) * steep / e50 * rho / (1 + rho)^2
  nested_log <- log(dose) - log(e50)
  d_steep <- (emax - e0) * nested_log * rho / (1 + rho)^2

  rbind(
    matrix(resid * d_lower, nrow = 1),
    matrix(resid * d_upper, nrow = 1),
    matrix(resid * d_ed50, nrow = 1),
    matrix(resid * d_steep, nrow = 1)
  )
}

#' Estimating equation for delta-effective dose (log-logistic)
#'
#' Computes the effective dose at level delta from the log-logistic model.
#' Should be stacked with [ee_loglogistic()].
#'
#' @param theta Numeric scalar. The ED(delta) parameter.
#' @param dose Numeric vector (used for dimension only).
#' @param delta Numeric effective dose level of interest.
#' @param lower Numeric lower limit parameter.
#' @param upper Numeric upper limit parameter.
#' @param ed50 Numeric estimated ED50.
#' @param steepness Numeric steepness parameter.
#'
#' @returns A 1-by-n matrix.
#'
#' @examples
#' # The lower limit is held at zero instead of being estimated: root length
#' # cannot fall below zero, and the five-parameter stack diverges on these
#' # data. The lower-limit row of the log-logistic equation is therefore
#' # dropped, leaving the upper limit, ED50, and steepness to be estimated
#' # alongside the ED90.
#' psi <- function(theta) {
#'   loglogistic <- ee_loglogistic(
#'     c(0, theta[1:3]),
#'     dose = inderjit$dose,
#'     response = inderjit$response
#'   )
#'   ed90 <- ee_loglogistic_ed(
#'     theta[4],
#'     dose = inderjit$dose,
#'     delta = 0.9,
#'     lower = 0,
#'     upper = theta[1],
#'     ed50 = theta[2],
#'     steepness = theta[3]
#'   )
#'   rbind(loglogistic[-1, , drop = FALSE], ed90)
#' }
#'
#' # This stacked system is sensitive to its starting values, so a reasonable
#' # init matters more here than it does for most equations.
#' m <- m_estimate(stacked_equations = psi, init = c(8, 2, 1, 5))
#'
#' # Upper limit, ED50, steepness, and the ED90
#' coef(m)
#'
#' @export
ee_loglogistic_ed <- function(
  theta,
  dose,
  delta,
  lower,
  upper,
  ed50,
  steepness
) {
  n <- length(dose)

  # Target response level depends on direction

  if (steepness < 0) {
    level <- upper * delta - lower * (1 - delta)
  } else {
    level <- upper * (1 - delta) - lower * delta
  }

  rho <- (theta[1] / ed50)^steepness
  pred <- lower + (upper - lower) / (1 + rho)

  matrix(rep(pred - level, n), nrow = 1)
}

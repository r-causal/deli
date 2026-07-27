#' Estimating equation for Rogan-Gladen correction
#'
#' Corrects for mismeasured binary outcomes using external validation data
#' to estimate sensitivity and specificity.
#'
#' @param theta Numeric vector of length 4: corrected proportion,
#'   naive proportion, sensitivity, specificity.
#' @param y Numeric vector of gold-standard measurements (only available
#'   in external validation sample where \code{r = 0}).
#' @param y_star Numeric vector of mismeasured outcome values (all observations).
#' @param r Numeric indicator: 1 for main study data, 0 for external validation.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A 4-by-n matrix.
#'
#' @examples
#' # A main study measures a binary outcome with an imperfect test. An external
#' # validation study measures both the test and the gold standard, and so
#' # informs the sensitivity and specificity used to correct the main study.
#' set.seed(2)
#' n_main <- 500
#' n_validation <- 400
#' n <- n_main + n_validation
#' y_true <- rbinom(n, 1, 0.25)
#' y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
#' r <- c(rep(1, n_main), rep(0, n_validation))
#'
#' # The gold standard is unobserved in the main study, so those positions carry
#' # a 0 placeholder. The placeholder never reaches an estimate: the sensitivity
#' # and specificity equations are multiplied by (1 - r).
#' y <- ifelse(r == 1, 0, y_true)
#'
#' psi <- function(theta) {
#'   ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
#' }
#'
#' m <- MEstimator(
#'   stacked_equations = psi,
#'   init = c(0.5, 0.5, 0.75, 0.75)
#' ) |>
#'   estimate()
#'
#' # Corrected prevalence, naive prevalence, sensitivity, specificity
#' coef(m)
#'
#' @export
ee_rogan_gladen <- function(theta, y, y_star, r, weights = NULL) {
  y <- as.numeric(y)
  y_star <- as.numeric(y_star)
  r <- as.numeric(r)
  n <- length(y_star)
  w <- generate_weights(n, weights)

  # Replace y with placeholder for main study obs
  y <- ifelse(r == 1, -999, y)

  mu <- theta[1] # Corrected proportion
  mu_star <- theta[2] # Naive proportion
  sens <- theta[3] # Sensitivity
  spec <- theta[4] # Specificity

  # Corrected mean: mu*(sens + spec - 1) = mu_star + spec - 1
  ef_corrected <- matrix(
    mu * (sens + spec - 1) - (mu_star + spec - 1),
    nrow = 1,
    ncol = n
  )

  # Naive mean: E[Y*] among main study (r=1)
  ef_naive <- matrix(r * w * (y_star - mu_star), nrow = 1)

  # Sensitivity: P(Y*=1 | Y=1) among validation (r=0)
  ef_sens <- matrix((y_star - sens) * (1 - r) * y * w, nrow = 1)

  # Specificity: P(Y*=0 | Y=0) among validation (r=0)
  ef_spec <- matrix((1 - y_star - spec) * (1 - r) * (1 - y) * w, nrow = 1)

  rbind(ef_corrected, ef_naive, ef_sens, ef_spec)
}

#' Estimating equation for extended Rogan-Gladen correction
#'
#' Extended version that conditions sensitivity and specificity on covariates
#' using logistic regression models.
#'
#' @param theta Numeric vector of length `1 + 2*p`: corrected proportion,
#'   then `p` sensitivity model parameters, then `p` specificity model
#'   parameters.
#' @param y Numeric vector of gold-standard measurements (validation sample
#'   where \code{r = 0}).
#' @param y_star Numeric vector of mismeasured outcome values.
#' @param r Numeric indicator: 1 for main study, 0 for validation.
#' @param X Numeric n-by-p design matrix for sensitivity/specificity models.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A `(1+2*p)`-by-n matrix.
#'
#' @examples
#' # The same validation design as ee_rogan_gladen(), with sensitivity and
#' # specificity now modeled by logistic regression. The design matrix here is
#' # intercept only, so both models estimate a single log-odds.
#' set.seed(123)
#' n <- 500
#' y_true <- rbinom(n, 1, 0.3)
#' y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
#' r <- c(rep(0, 200), rep(1, 300))
#'
#' # Gold standard observed only in the validation sample (r == 0)
#' y <- ifelse(r == 0, y_true, 0)
#' X <- cbind(rep(1, n))
#'
#' psi <- function(theta) {
#'   ee_rogan_gladen_extended(theta, y = y, y_star = y_star, r = r, X = X)
#' }
#'
#' m <- MEstimator(stacked_equations = psi, init = c(0.5, 1, 1)) |>
#'   estimate()
#'
#' # Corrected prevalence, then the sensitivity and specificity intercepts
#' coef(m)
#'
#' @export
ee_rogan_gladen_extended <- function(theta, y, y_star, r, X, weights = NULL) {
  y <- as.numeric(y)
  y_star <- as.numeric(y_star)
  r <- as.numeric(r)
  X <- as.matrix(X)
  n <- length(y_star)
  p <- ncol(X)
  w <- generate_weights(n, weights)

  # Replace y with placeholder for main study obs
  y <- ifelse(r == 1, -999, y)

  mu <- theta[1]
  sens_params <- theta[2:(1 + p)]
  spec_params <- theta[(2 + p):(1 + 2 * p)]

  # Sensitivity model: logistic regression of y_star on X among Y=1, r=0
  ee_sens <- ee_regression(
    theta = sens_params,
    X = X,
    y = y_star,
    model = "logistic",
    weights = weights
  ) *
    matrix((1 - r) * y, nrow = p, ncol = n, byrow = TRUE)

  sens_i <- inverse_logit(pt_as_vector(X %*% sens_params))

  # Specificity model: logistic regression of (1-y_star) on X among Y=0, r=0
  ee_spec <- ee_regression(
    theta = spec_params,
    X = X,
    y = 1 - y_star,
    model = "logistic",
    weights = weights
  ) *
    matrix((1 - r) * (1 - y), nrow = p, ncol = n, byrow = TRUE)

  spec_i <- inverse_logit(pt_as_vector(X %*% spec_params))

  # Individual-level Rogan-Gladen correction
  rg_equation <- (y_star + spec_i - 1) / (sens_i + spec_i - 1)
  ef_corrected <- matrix(r * (rg_equation - mu) * w, nrow = 1)

  rbind(ef_corrected, ee_sens, ee_spec)
}

#' Estimating equation for regression calibration
#'
#' Corrects for measurement error in a binary predictor using external
#' validation data. Scales the naive coefficient by the calibration factor.
#'
#' @param theta Numeric vector: the corrected coefficient first, then the
#'   calibration coefficients for the design `cbind(a_star, X)` (the `a_star`
#'   coefficient first). Length `2 + ncol(X)` when `X` is supplied, or length
#'   `3` when `X = NULL` (an intercept-only calibration model).
#' @param beta Numeric scalar. External estimate of the coefficient for the
#'   mismeasured predictor on the outcome.
#' @param a Numeric vector of gold-standard action measurements (validation
#'   sample only, where \code{r = 0}).
#' @param a_star Numeric vector of mismeasured action values.
#' @param r Numeric indicator: 0 for validation, 1 for main study.
#' @param X Optional design matrix for calibration model. Default `NULL`
#'   uses intercept only.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A `length(theta)`-by-n matrix (`2 + ncol(X)` rows, or `3` rows
#'   when `X = NULL`).
#'
#' @examples
#' # A binary exposure is measured with error in the main study. The external
#' # validation study regresses the gold-standard exposure on the mismeasured
#' # one, and that calibration slope rescales the naive outcome coefficient.
#' set.seed(789)
#' n <- 500
#' a_true <- rbinom(n, 1, 0.5)
#' a_star <- ifelse(a_true == 1, rbinom(n, 1, 0.85), rbinom(n, 1, 0.1))
#' r <- c(rep(0, 200), rep(1, 300))
#'
#' # Gold standard observed only in the validation sample (r == 0)
#' a <- ifelse(r == 0, a_true, 0)
#'
#' # `beta` is the naive coefficient for the mismeasured exposure, supplied here
#' # as a fixed external value. Stack an outcome model and pass its coefficient
#' # instead to propagate the uncertainty in that estimate as well.
#' psi <- function(theta) {
#'   ee_regression_calibration(theta, beta = 0.8, a = a, a_star = a_star, r = r)
#' }
#'
#' m <- MEstimator(stacked_equations = psi, init = c(1, 0.1, 0.5)) |>
#'   estimate()
#'
#' # Corrected coefficient, then the calibration slope and intercept
#' coef(m)
#'
#' @export
ee_regression_calibration <- function(
  theta,
  beta,
  a,
  a_star,
  r,
  X = NULL,
  weights = NULL
) {
  a <- as.numeric(a)
  a_star <- as.numeric(a_star)
  r <- as.numeric(r)
  n <- length(a_star)
  w <- generate_weights(n, weights)

  # Replace a with placeholder for main study obs
  a <- ifelse(r == 1, -999, a)

  # Default to a_star + intercept design matrix (matches Python column order)
  if (is.null(X)) {
    Xcal <- cbind(a_star, 1)
  } else {
    Xcal <- cbind(a_star, as.matrix(X))
  }
  p <- ncol(Xcal)

  beta_corrected <- theta[1]
  gamma <- theta[2:(1 + p)]

  # Corrected coefficient: beta / gamma[1] (a_star coefficient, first column)
  ef_corrected <- matrix(
    beta / gamma[1] - beta_corrected,
    nrow = 1,
    ncol = n
  )

  # Calibration model: linear model of a ~ a_star (+ X) using validation data
  eta <- pt_as_vector(Xcal %*% gamma)
  resid <- a - eta
  ef_cal <- t(Xcal * (w * (1 - r) * resid))

  rbind(ef_corrected, ef_cal)
}

#' Standard S3 generics for deli estimators
#'
#' Methods for base R generics [stats::coef()], [stats::vcov()],
#' [stats::confint()], and [stats::nobs()] so that deli estimator objects
#' interoperate with the broader R modeling ecosystem.
#'
#' @param object A fitted `MEstimator` or `GMMEstimator` object (after calling
#'   [estimate()]).
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level The confidence level required. Default `0.95`.
#' @param ... Additional arguments (currently unused).
#'
#' @returns
#' - `coef()`: Named numeric vector of parameter estimates.
#' - `vcov()`: Named variance-covariance matrix.
#' - `confint()`: Matrix with columns `"lower"` and `"upper"`.
#' - `nobs()`: Integer number of observations.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' coef(fit)
#'
#' vcov(fit)
#'
#' confint(fit)
#'
#' nobs(fit)
#'
#' @name deli-generics
#' @importFrom stats coef vcov confint nobs
NULL

# ---- External generic declarations ------------------------------------------

stats_coef <- new_external_generic("stats", "coef", "object")
stats_vcov <- new_external_generic("stats", "vcov", "object")
stats_confint <- new_external_generic("stats", "confint", "object")
stats_nobs <- new_external_generic("stats", "nobs", "object")

# ---- coef --------------------------------------------------------------------

method(stats_coef, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@theta
}

# ---- vcov --------------------------------------------------------------------

method(stats_vcov, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@variance
}

# ---- confint -----------------------------------------------------------------

method(stats_confint, deli_estimator) <- function(
  object,
  parm,
  level = 0.95,
  ...
) {
  check_estimated(object)
  ci <- confidence_intervals(object, alpha = 1 - level)
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}

# ---- nobs --------------------------------------------------------------------

method(stats_nobs, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@n_obs
}

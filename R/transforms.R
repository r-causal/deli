#' Logistic transformation
#'
#' Transforms probabilities into log-odds: \eqn{\log(p / (1 - p))}.
#'
#' @param prob A numeric value or vector of probabilities.
#'
#' @return Numeric log-odds values.
#'
#' @export
#' @examples
#' logit(0.5)
#' logit(c(0.1, 0.5, 0.9))
logit <- function(prob) {
  # Written as log(prob / (1 - prob)) rather than qlogis() so the arithmetic and
  # Math group generics carry tangents through automatically under exact-mode
  # autodiff, mirroring inverse_logit() and the Python implementation. Numerics
  # are identical to qlogis() for numeric input.
  log(prob / (1 - prob))
}

#' Inverse logistic transformation
#'
#' Transforms log-odds into probabilities: \eqn{1 / (1 + \exp(-x))}.
#'
#' @param logodds A numeric value or vector of log-odds.
#'
#' @return Numeric probability values.
#'
#' @export
#' @examples
#' inverse_logit(0)
#' inverse_logit(c(-2, 0, 2))
inverse_logit <- function(logodds) {
  # Written as 1 / (1 + exp(-x)) rather than plogis() so the arithmetic and
  # Math group generics carry tangents through automatically, mirroring the
  # Python implementation in delicatessen/utilities.py.
  1 / (1 + exp(-logodds))
}

#' Identity transformation
#'
#' Returns the input unchanged. Used as a no-op transformation in contexts
#' that accept an arbitrary transformation function.
#'
#' @param value Any value.
#'
#' @return The input `value`, unchanged.
#'
#' @export
#' @examples
#' identity_transform(42)
#' identity_transform(c(1, 2, 3))
identity_transform <- function(value) {
  value
}

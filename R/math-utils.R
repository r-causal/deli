#' Polygamma function
#'
#' Wrapper around [psigamma()] with argument order matching Python's
#' `scipy.special.polygamma(n, x)`.
#'
#' @param n Integer order of the derivative of the digamma function.
#' @param x Numeric value or vector.
#'
#' @return Numeric polygamma values.
#'
#' @export
#' @examples
#' deli_polygamma(0, 1)
#' deli_polygamma(1, c(1, 2, 5))
deli_polygamma <- function(n, x) {
  # A whole parameter vector arrives as a PrimalTangentVector under exact mode;
  # normalize it to a tangent array so the array branch below differentiates it
  # elementwise rather than falling through to the numeric kernel.
  if (inherits(x, "PrimalTangentVector")) {
    parts <- pt_arrays(x)
    x <- primal_tangent_array(parts$primal, parts$tangent)
  }
  # Dispatch on tangent-carrying inputs before the numeric kernel: the tangent
  # of polygamma(n, x) is polygamma(n + 1, x), matching Python's PrimalTangentPairs
  if (is_pt(x)) {
    primal <- psigamma(x$primal, deriv = n)
    tangent <- x$tangent * psigamma(x$primal, deriv = n + 1)
    return(primal_tangent(primal, tangent))
  }
  if (is_pt_array(x)) {
    primal <- psigamma(x$primal, deriv = n)
    tangent <- x$tangent * psigamma(x$primal, deriv = n + 1)
    return(primal_tangent_array(primal, tangent))
  }
  psigamma(x, deriv = n)
}

#' Digamma function
#'
#' Wrapper around [base::digamma()].
#'
#' @param z Numeric value or vector.
#'
#' @return Numeric digamma values.
#'
#' @export
#' @examples
#' deli_digamma(1)
#' deli_digamma(c(0.5, 1, 2))
deli_digamma <- function(z) {
  # A whole parameter vector arrives as a PrimalTangentVector under exact mode;
  # normalize it to a tangent array so the array branch below differentiates it
  # elementwise rather than falling through to the numeric kernel.
  if (inherits(z, "PrimalTangentVector")) {
    parts <- pt_arrays(z)
    z <- primal_tangent_array(parts$primal, parts$tangent)
  }
  # Dispatch on tangent-carrying inputs before the as.numeric() coercion below,
  # which would otherwise strip the tangent. The tangent of digamma(z) is
  # trigamma(z), matching Python's PrimalTangentPairs.
  if (is_pt(z)) {
    return(primal_tangent(digamma(z$primal), z$tangent * trigamma(z$primal)))
  }
  if (is_pt_array(z)) {
    return(primal_tangent_array(
      digamma(z$primal),
      z$tangent * trigamma(z$primal)
    ))
  }
  # Return NaN directly for non-positive integers (poles of digamma)
  # without triggering base R's NaN warning
  z <- as.numeric(z)
  out <- numeric(length(z))
  poles <- z <= 0 & z == trunc(z)
  out[poles] <- NaN
  if (any(!poles)) {
    out[!poles] <- digamma(z[!poles])
  }
  out
}

#' Standard normal CDF
#'
#' Evaluates the cumulative distribution function of the standard normal
#' distribution. Wrapper around [pnorm()].
#'
#' @param x Numeric value or vector of quantiles.
#'
#' @return Numeric CDF values.
#'
#' @export
#' @examples
#' standard_normal_cdf(0)
#' standard_normal_cdf(c(-1.96, 0, 1.96))
standard_normal_cdf <- function(x) {
  # A whole parameter vector arrives as a PrimalTangentVector under exact mode;
  # normalize it to a tangent array so the array branch below differentiates it
  # elementwise rather than falling through to the numeric kernel.
  if (inherits(x, "PrimalTangentVector")) {
    parts <- pt_arrays(x)
    x <- primal_tangent_array(parts$primal, parts$tangent)
  }
  # Dispatch on tangent-carrying inputs: the tangent of the CDF is the PDF,
  # matching Python's PrimalTangentPairs.normal_cdf
  if (is_pt(x)) {
    return(primal_tangent(pnorm(x$primal), x$tangent * dnorm(x$primal)))
  }
  if (is_pt_array(x)) {
    return(primal_tangent_array(pnorm(x$primal), x$tangent * dnorm(x$primal)))
  }
  pnorm(x)
}

#' Standard normal PDF
#'
#' Evaluates the probability density function of the standard normal
#' distribution. Wrapper around [dnorm()].
#'
#' @param x Numeric value or vector of quantiles.
#'
#' @return Numeric density values.
#'
#' @export
#' @examples
#' standard_normal_pdf(0)
#' standard_normal_pdf(c(-1, 0, 1))
standard_normal_pdf <- function(x) {
  # A whole parameter vector arrives as a PrimalTangentVector under exact mode;
  # normalize it to a tangent array so the array branch below differentiates it
  # elementwise rather than falling through to the numeric kernel.
  if (inherits(x, "PrimalTangentVector")) {
    parts <- pt_arrays(x)
    x <- primal_tangent_array(parts$primal, parts$tangent)
  }
  # Dispatch on tangent-carrying inputs: the tangent of the PDF is -x times the
  # PDF, matching Python's PrimalTangentPairs.normal_pdf
  if (is_pt(x)) {
    return(primal_tangent(
      dnorm(x$primal),
      x$tangent * (-x$primal * dnorm(x$primal))
    ))
  }
  if (is_pt_array(x)) {
    return(primal_tangent_array(
      dnorm(x$primal),
      x$tangent * (-x$primal * dnorm(x$primal))
    ))
  }
  dnorm(x)
}

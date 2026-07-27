#' Polygamma function
#'
#' @description
#' Evaluates the `n`th derivative of the digamma function.
#'
#' This is the equivalent of [base::psigamma()] and returns identical values for
#' numeric input, except that it carries derivatives. Exact differentiation
#' (`deriv_method = "exact"`) propagates a tangent alongside each value, and
#' `deli_polygamma()` recognizes a tangent-carrying argument and applies the
#' analytic rule itself, the polygamma function of order `n + 1`. `psigamma()`
#' hands its argument to compiled code without dispatching, and errors on such
#' an argument. Use `deli_polygamma()` inside estimating equations and inside
#' transforms passed to [delta_method()], and `psigamma()` for ordinary numeric
#' work.
#'
#' **The two take their arguments in opposite orders.** `deli_polygamma(n, x)`
#' takes the derivative order first, matching Python's
#' `scipy.special.polygamma(n, x)`; `psigamma(x, deriv = n)` takes it second.
#' `deli_polygamma(1, 3)` is `psigamma(3, deriv = 1)`, not
#' `psigamma(1, deriv = 3)`. A positional substitution between the two computes
#' a different quantity and raises no error.
#'
#' @param n Integer order of the derivative of the digamma function.
#' @param x Numeric value or vector.
#'
#' @returns Numeric polygamma values.
#'
#' @inheritSection inverse_logit Exact differentiation
#'
#' @seealso [base::psigamma()], the base R equivalent for ordinary numeric work,
#'   which takes its arguments in the other order, and [deli_digamma()] for the
#'   order-zero case.
#'
#' @examples
#' deli_polygamma(0, 1)
#' deli_polygamma(1, c(1, 2, 5))
#'
#' # The same values as the base R counterpart, with the arguments swapped
#' all.equal(deli_polygamma(1, c(1, 2, 5)), psigamma(c(1, 2, 5), deriv = 1))
#'
#' # A gamma GLM, whose last parameter is the log of the shape
#' m <- m_estimate(
#'   mpg ~ wt,
#'   data = mtcars,
#'   .ee = ee_glm,
#'   distribution = "gamma",
#'   link = "log",
#'   init = c(3, 0, 0)
#' )
#'
#' # Variance of the trigamma function at the estimated shape. Writing
#' # `psigamma(exp(theta[3]), deriv = 1)` here instead would error, because
#' # exact differentiation hands the transform a tangent-carrying argument.
#' delta_method(
#'   m,
#'   transform = function(theta) deli_polygamma(1, exp(theta[3])),
#'   deriv_method = "exact"
#' )
#'
#' @export
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
#' @description
#' Evaluates the digamma function, the first derivative of the log gamma
#' function.
#'
#' This is the equivalent of [base::digamma()], and unlike most of the deli
#' utilities it is not required anywhere: `digamma()` belongs to the `Math`
#' group generic and has a tangent rule, so it differentiates under
#' `deriv_method = "exact"` just as `deli_digamma()` does. The two agree value
#' for value, including at the poles and on missing input. The only difference
#' is that `deli_digamma()` returns `NaN` at the non-positive integers, where
#' the digamma function has poles, without base R's `NaNs produced` warning.
#' `deli_digamma()` exists so that code translated from Python delicatessen,
#' where the function is called `digamma()`, can keep its shape.
#'
#' @param z Numeric value or vector.
#'
#' @returns Numeric digamma values.
#'
#' @seealso [base::digamma()], which is usable in the same positions, and
#'   [deli_polygamma()] for the higher-order derivatives. [inverse_logit()] has
#'   the full list of deli utilities and their base R counterparts, including
#'   the five that base R cannot stand in for under exact differentiation.
#'
#' @examples
#' deli_digamma(1)
#' deli_digamma(c(0.5, 1, 2))
#'
#' # The same values as the base R counterpart
#' all.equal(deli_digamma(c(0.5, 1, 2)), digamma(c(0.5, 1, 2)))
#'
#' # NaN at the poles, where base R warns, and missing values pass through
#' deli_digamma(c(-1, 0, 2.5, NA))
#'
#' @export
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
  # Return NaN directly for non-positive integers (poles of digamma) without
  # triggering base R's NaN warning. A missing input must not reach the pole
  # test: `NA <= 0` is NA, which would poison both the logical index and the
  # subscripted assignment. Guard the poles against NA/NaN, propagate the missing
  # values unchanged (preserving the NA vs NaN distinction, matching
  # base::digamma), and evaluate digamma only on the finite non-pole entries.
  z <- as.numeric(z)
  out <- numeric(length(z))
  na_z <- is.na(z)
  poles <- !na_z & z <= 0 & z == trunc(z)
  out[poles] <- NaN
  out[na_z] <- z[na_z]
  finite_regular <- !poles & !na_z
  out[finite_regular] <- digamma(z[finite_regular])
  out
}

#' Standard normal CDF
#'
#' @description
#' Evaluates the cumulative distribution function of the standard normal
#' distribution.
#'
#' This is the equivalent of [stats::pnorm()] at the default mean and standard
#' deviation, and returns identical values for numeric input, except that it
#' carries derivatives. Exact differentiation (`deriv_method = "exact"`)
#' propagates a tangent alongside each value, and `standard_normal_cdf()`
#' recognizes a tangent-carrying argument and applies the analytic rule itself,
#' the standard normal density. `pnorm()` hands its argument to compiled code
#' without dispatching, and errors on such an argument. Use
#' `standard_normal_cdf()` inside estimating equations and inside transforms
#' passed to [delta_method()], and `pnorm()` for ordinary numeric work.
#'
#' @param x Numeric value or vector of quantiles.
#'
#' @returns Numeric CDF values.
#'
#' @inheritSection inverse_logit Exact differentiation
#'
#' @seealso [stats::pnorm()], the base R equivalent for ordinary numeric work,
#'   and [standard_normal_pdf()] for the density.
#'
#' @examples
#' standard_normal_cdf(0)
#' standard_normal_cdf(c(-1.96, 0, 1.96))
#'
#' # The same values as the base R counterpart
#' all.equal(standard_normal_cdf(c(-1.96, 0, 1.96)), pnorm(c(-1.96, 0, 1.96)))
#'
#' # A probit regression, whose mean model is the standard normal CDF
#' m <- m_estimate(
#'   vs ~ mpg,
#'   data = mtcars,
#'   .ee = ee_glm,
#'   distribution = "binomial",
#'   link = "probit"
#' )
#'
#' # Variance of the fitted probability at mpg = 20. Writing
#' # `pnorm(theta[1] + theta[2] * 20)` here instead would error, because exact
#' # differentiation hands the transform a tangent-carrying argument.
#' delta_method(
#'   m,
#'   transform = function(theta) standard_normal_cdf(theta[1] + theta[2] * 20),
#'   deriv_method = "exact"
#' )
#'
#' @export
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
#' @description
#' Evaluates the probability density function of the standard normal
#' distribution.
#'
#' This is the equivalent of [stats::dnorm()] at the default mean and standard
#' deviation, and returns identical values for numeric input, except that it
#' carries derivatives. Exact differentiation (`deriv_method = "exact"`)
#' propagates a tangent alongside each value, and `standard_normal_pdf()`
#' recognizes a tangent-carrying argument and applies the analytic rule itself,
#' \eqn{-x} times the density. `dnorm()` hands its argument to compiled code
#' without dispatching, and errors on such an argument. Use
#' `standard_normal_pdf()` inside estimating equations and inside transforms
#' passed to [delta_method()], and `dnorm()` for ordinary numeric work.
#'
#' @param x Numeric value or vector of quantiles.
#'
#' @returns Numeric density values.
#'
#' @inheritSection inverse_logit Exact differentiation
#'
#' @seealso [stats::dnorm()], the base R equivalent for ordinary numeric work,
#'   and [standard_normal_cdf()] for the distribution function.
#'
#' @examples
#' standard_normal_pdf(0)
#' standard_normal_pdf(c(-1, 0, 1))
#'
#' # The same values as the base R counterpart
#' all.equal(standard_normal_pdf(c(-1, 0, 1)), dnorm(c(-1, 0, 1)))
#'
#' m <- m_estimate(
#'   vs ~ mpg,
#'   data = mtcars,
#'   .ee = ee_glm,
#'   distribution = "binomial",
#'   link = "probit"
#' )
#'
#' # Variance of the marginal effect of mpg at mpg = 20, the standard normal
#' # density at the linear predictor times the coefficient. Writing
#' # `dnorm(theta[1] + theta[2] * 20)` here instead would error, because exact
#' # differentiation hands the transform a tangent-carrying argument.
#' delta_method(
#'   m,
#'   transform = function(theta) {
#'     standard_normal_pdf(theta[1] + theta[2] * 20) * theta[2]
#'   },
#'   deriv_method = "exact"
#' )
#'
#' @export
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

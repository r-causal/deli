#' Logistic transformation
#'
#' @description
#' Transforms probabilities into log-odds: \eqn{\log(p / (1 - p))}.
#'
#' This is the equivalent of [stats::qlogis()] and returns identical values for
#' numeric input, except that it carries derivatives. Exact differentiation
#' (`deriv_method = "exact"`) propagates a tangent alongside each value through
#' the arithmetic and `Math` group generics. `logit()` is written as
#' `log(prob / (1 - prob))`, so those group generic methods carry a tangent
#' through it; `qlogis()` hands its argument to compiled code without
#' dispatching, and errors on a tangent-carrying argument. Use `logit()`
#' inside estimating equations and inside transforms passed to [delta_method()],
#' and `qlogis()` for ordinary numeric work.
#'
#' @param prob A numeric value or vector of probabilities.
#'
#' @returns Numeric log-odds values.
#'
#' @inheritSection inverse_logit Exact differentiation
#'
#' @seealso [stats::qlogis()], the base R equivalent for ordinary numeric work,
#'   and [inverse_logit()], which inverts this transformation.
#'
#' @examples
#' logit(0.5)
#' logit(c(0.1, 0.5, 0.9))
#'
#' # The same values as the base R counterpart
#' all.equal(logit(c(0.1, 0.5, 0.9)), qlogis(c(0.1, 0.5, 0.9)))
#'
#' # The proportion of cars with a manual transmission
#' m <- m_estimate(
#'   stacked_equations = function(theta) ee_mean(theta, y = mtcars$am),
#'   init = 0.5
#' )
#'
#' # Variance of the log-odds of that proportion. Writing `qlogis(theta[1])`
#' # here instead would error, because exact differentiation hands the
#' # transform a tangent-carrying argument.
#' delta_method(
#'   m,
#'   transform = function(theta) logit(theta[1]),
#'   deriv_method = "exact"
#' )
#'
#' @export
logit <- function(prob) {
  # Written as log(prob / (1 - prob)) rather than qlogis() so the arithmetic and
  # Math group generics carry tangents through automatically under exact-mode
  # autodiff, mirroring inverse_logit() and the Python implementation. Numerics
  # are identical to qlogis() for numeric input.
  log(prob / (1 - prob))
}

#' Inverse logistic transformation
#'
#' @description
#' Transforms log-odds into probabilities: \eqn{1 / (1 + \exp(-x))}.
#'
#' This is the equivalent of [stats::plogis()] and returns identical values for
#' numeric input, except that it carries derivatives. Exact differentiation
#' (`deriv_method = "exact"`) propagates a tangent alongside each value through
#' the arithmetic and `Math` group generics. `inverse_logit()` is written as
#' `1 / (1 + exp(-logodds))`, so those group generic methods carry a tangent
#' through it; `plogis()` hands its argument to compiled code without
#' dispatching, and errors on a tangent-carrying argument. Use
#' `inverse_logit()` inside estimating equations and inside transforms passed to
#' [delta_method()], and `plogis()` for ordinary numeric work.
#'
#' @param logodds A numeric value or vector of log-odds.
#'
#' @returns Numeric probability values.
#'
#' @section Exact differentiation:
#' `deriv_method = "exact"` is forward-mode automatic differentiation: it
#' replaces each value with an object carrying both the value and its
#' derivative. deli registers `Ops`, `Math`, and `Summary` group generic methods
#' for those objects, and a function differentiates exactly only if it
#' dispatches to one of them. `plogis()`, `qlogis()`, `pnorm()`, `dnorm()`, and
#' `psigamma()` never dispatch: each hands its argument straight to compiled
#' code through `.Call()` or `.Internal()`, which requires a plain number, so
#' the call stops with `Non-numeric argument to mathematical function`.
#'
#' Group membership is necessary but not sufficient. deli's `Math` method
#' implements tangent rules for part of the group and raises an error for the
#' rest: `gamma()`, `round()`, `signif()`, `trunc()`, `cumprod()`, `cummax()`,
#' and `cummin()`. `vignette("getting-started")` lists the whole supported
#' surface.
#'
#' Each deli utility below returns the same values as its base R counterpart for
#' numeric input. What separates them is whether the counterpart survives exact
#' mode:
#'
#' | deli function | base R counterpart | base R under `deriv_method = "exact"` |
#' | --- | --- | --- |
#' | [inverse_logit()] | [stats::plogis()] | errors |
#' | [logit()] | [stats::qlogis()] | errors |
#' | [standard_normal_cdf()] | [stats::pnorm()] | errors |
#' | [standard_normal_pdf()] | [stats::dnorm()] | errors |
#' | [deli_polygamma()] | [base::psigamma()] | errors |
#' | [deli_digamma()] | [base::digamma()] | works |
#' | [identity_transform()] | [base::identity()] | works |
#'
#' For the first five rows, use the deli function inside estimating equations
#' and inside transforms passed to [delta_method()], and the base R function
#' everywhere else: simulating data, post-fit display, plain numeric work. For
#' the last two rows the base R function is usable everywhere, because
#' `digamma()` is a `Math` group member with a tangent rule and `identity()`
#' passes its argument through untouched.
#'
#' The polygamma row is the one place where the arguments do not line up.
#' `deli_polygamma(n, x)` takes the derivative order first and
#' `psigamma(x, deriv = n)` takes it second, so a positional substitution
#' between the two computes a different quantity and raises no error.
#'
#' @seealso [stats::plogis()], the base R equivalent for ordinary numeric work,
#'   and [logit()], which inverts this transformation.
#'
#' @examples
#' inverse_logit(0)
#' inverse_logit(c(-2, 0, 2))
#'
#' # The same values as the base R counterpart
#' all.equal(inverse_logit(c(-2, 0, 2)), plogis(c(-2, 0, 2)))
#'
#' m <- m_estimate(
#'   vs ~ mpg,
#'   data = mtcars,
#'   .ee = ee_regression,
#'   model = "logistic"
#' )
#'
#' # Variance of the fitted probability at mpg = 20. Writing
#' # `plogis(theta[1] + theta[2] * 20)` here instead would error, because exact
#' # differentiation hands the transform a tangent-carrying argument.
#' delta_method(
#'   m,
#'   transform = function(theta) inverse_logit(theta[1] + theta[2] * 20),
#'   deriv_method = "exact"
#' )
#'
#' @export
inverse_logit <- function(logodds) {
  # Written as 1 / (1 + exp(-x)) rather than plogis() so the arithmetic and
  # Math group generics carry tangents through automatically, mirroring the
  # Python implementation in delicatessen/utilities.py.
  1 / (1 + exp(-logodds))
}

#' Identity transformation
#'
#' @description
#' Returns the input unchanged. Used as a no-op transformation in contexts
#' that accept an arbitrary transformation function.
#'
#' This is [base::identity()] under a name that does not mask it. The two are
#' interchangeable everywhere, including under exact differentiation
#' (`deriv_method = "exact"`), since returning the argument untouched also
#' returns any tangent it carries. `identity_transform()` exists so that code
#' translated from Python delicatessen, where the function is called
#' `identity()`, can keep its shape. In R code, prefer `identity()`.
#'
#' @param value Any value.
#'
#' @returns The input `value`, unchanged.
#'
#' @seealso [base::identity()], which does the same thing and is the idiomatic
#'   choice in R. [inverse_logit()] has the full list of deli utilities and
#'   their base R counterparts.
#'
#' @examples
#' identity_transform(42)
#' identity_transform(c(1, 2, 3))
#'
#' # The same value as base R, and the same object
#' identical(identity_transform(c(1, 2, 3)), identity(c(1, 2, 3)))
#'
#' @export
identity_transform <- function(value) {
  value
}

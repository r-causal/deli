#' Confidence intervals for M-Estimator parameters
#'
#' @description
#' Computes two-sided Wald-type \eqn{(1 - \alpha) \times 100\%} confidence
#' intervals using the point estimates and sandwich variance:
#' \eqn{\hat{\theta} \pm c_{\alpha/2} \times \widehat{SE}(\hat{\theta})}.
#' The critical value \eqn{c_{\alpha/2}} comes from the standard normal
#' distribution by default. When a `finite_correction` is set on the fit, it
#' comes instead from the t-distribution with \eqn{n - p} degrees of freedom,
#' matching the finite-sample adjustment of the variance.
#'
#' This function mirrors `m.confidence_intervals()` in Python delicatessen, so
#' code translated from Python can keep its shape.
#'
#' @param object A fitted `MEstimator` object (after calling [estimate()]).
#' @param alpha Numeric significance level, between 0 and 1. Default `0.05`
#'   for 95% confidence intervals.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A p-by-2 matrix with columns `"lower"` and `"upper"`.
#'
#' @seealso [`confint()`][deli-generics], the standard R accessor for the same
#'   intervals. It returns identical values but is parameterized by
#'   `level = 0.95` where this function takes `alpha = 0.05`.
#'
#' @export
#' @examples
#' psi <- function(theta) {
#'   y <- c(1, 2, 3, 4, 5)
#'   matrix(y - theta[1], nrow = 1)
#' }
#' m <- m_estimate(stacked_equations = psi, init = 0)
#' confidence_intervals(m)
confidence_intervals <- new_generic(
  "confidence_intervals",
  "object",
  function(object, alpha = 0.05, ...) {
    S7::S7_dispatch()
  }
)

method(confidence_intervals, deli_estimator) <- function(
  object,
  alpha = 0.05,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)
  check_alpha(alpha)

  # Critical value from normal, or t under a finite-sample correction
  cv <- critical_value(object, alpha)

  # Standard errors from diagonal of variance matrix
  se <- sqrt(diag(object@variance))

  # Wald-type confidence intervals
  lower <- object@theta - cv * se
  upper <- object@theta + cv * se

  ci <- cbind(lower = lower, upper = upper)
  rownames(ci) <- names(object@theta)
  ci
}

#' Z-scores for M-Estimator parameters
#'
#' @description
#' Computes Wald-type Z-scores:
#' \eqn{(\hat{\theta} - \theta_0) / \widehat{SE}(\hat{\theta})}.
#'
#' This function mirrors `m.z_scores()` in Python delicatessen, so code
#' translated from Python can keep its shape.
#'
#' @param object A fitted `MEstimator` object (after calling [estimate()]).
#' @param null Numeric null hypothesis value(s). Default `0`.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A numeric vector of Z-scores.
#'
#' @seealso [`summary()`][deli-display] and [`tidy()`][deli-tidiers], which
#'   report the same Z-scores in table form alongside the other parameter-level
#'   results.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' z_scores(fit)
#'
#' @export
z_scores <- new_generic("z_scores", "object", function(object, null = 0, ...) {
  S7::S7_dispatch()
})

method(z_scores, deli_estimator) <- function(object, null = 0, ...) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)
  check_null_length(null, length(object@theta))

  se <- sqrt(diag(object@variance))
  (object@theta - null) / se
}

#' P-values for M-Estimator parameters
#'
#' @description
#' Computes two-sided Wald-type P-values from Z-scores. The Z-scores are
#' compared to the standard normal distribution by default. When a
#' `finite_correction` is set on the fit, they are compared instead to the
#' t-distribution with \eqn{n - p} degrees of freedom.
#'
#' This function mirrors `m.p_values()` in Python delicatessen, so code
#' translated from Python can keep its shape.
#'
#' @param object A fitted `MEstimator` object (after calling [estimate()]).
#' @param null Numeric null hypothesis value(s). Default `0`.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A numeric vector of P-values.
#'
#' @seealso [`summary()`][deli-display] and [`tidy()`][deli-tidiers], which
#'   report the same P-values in table form alongside the other parameter-level
#'   results.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' p_values(fit)
#'
#' @export
p_values <- new_generic("p_values", "object", function(object, null = 0, ...) {
  S7::S7_dispatch()
})

method(p_values, deli_estimator) <- function(object, null = 0, ...) {
  rlang::check_dots_empty(call = rlang::caller_env())
  z <- z_scores(object, null = null)
  tail_probability(object, z)
}

#' S-values (surprisal) for M-Estimator parameters
#'
#' @description
#' Computes Shannon Information values (S-values) as \eqn{S = -\log_2(P)},
#' where \eqn{P} is the corresponding P-value.
#'
#' This function mirrors `m.s_values()` in Python delicatessen, so code
#' translated from Python can keep its shape.
#'
#' @details
#' A P-value small enough to underflow to exactly zero has an S-value of `Inf`.
#' That is the limit the surprisal is heading toward rather than a defect. The
#' smallest P a double can hold is \eqn{2^{-1074}}, so a P that arrives as
#' zero stands for more than a thousand bits of surprisal, past the range a
#' double can name. The infinity reports evidence beyond measurement, where any
#' finite substitute would name a number the fit does not support.
#'
#' @param object A fitted `MEstimator` object (after calling [estimate()]).
#' @param null Numeric null hypothesis value(s). Default `0`.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A numeric vector of S-values.
#'
#' @seealso [`summary()`][deli-display] and [`tidy()`][deli-tidiers], which
#'   report the same S-values in table form alongside the other parameter-level
#'   results.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' s_values(fit)
#'
#' @export
s_values <- new_generic("s_values", "object", function(object, null = 0, ...) {
  S7::S7_dispatch()
})

method(s_values, deli_estimator) <- function(object, null = 0, ...) {
  rlang::check_dots_empty(call = rlang::caller_env())
  p <- p_values(object, null = null)
  -log2(p)
}

# ---- Reference-distribution helpers ------------------------------------------

#' Two-sided critical value for confidence intervals
#'
#' Returns the normal quantile by default, or the t-distribution quantile with
#' `n_obs - n_params` degrees of freedom when the fit carries a finite-sample
#' correction. Mirrors delicatessen's `confidence_intervals` behavior.
#'
#' @noRd
critical_value <- function(object, alpha) {
  if (is.null(object@finite_correction)) {
    return(qnorm(1 - alpha / 2))
  }

  # Degrees of freedom for the t reference distribution
  df <- object@n_obs - object@n_params
  qt(1 - alpha / 2, df = df)
}

#' Two-sided tail probability for P-values
#'
#' Returns twice the upper-tail probability of the reference distribution: the
#' standard normal by default, or the t-distribution with `n_obs - n_params`
#' degrees of freedom under a finite-sample correction.
#'
#' @noRd
tail_probability <- function(object, z) {
  if (is.null(object@finite_correction)) {
    return(2 * pnorm(-abs(z)))
  }

  # Degrees of freedom for the t reference distribution
  df <- object@n_obs - object@n_params
  2 * pt(-abs(z), df = df)
}

#' Refuse an estimator that was never solved
#'
#' The requirement of everything that reports on the solve rather than on the
#' sandwich: the estimates, the counts taken with them, and the point
#' predictions formed from them. Each is recorded by the solve whatever became
#' of the variance, so each answers wherever there are estimates at all.
#'
#' @param object The estimator the caller asked something of.
#' @returns Invisible `NULL` where the fit carries estimates; otherwise it
#'   throws.
#' @noRd
check_has_estimates <- function(object) {
  if (is.null(object@theta)) {
    cli::cli_abort(
      "Cannot compute inference before calling {.fn estimate}.",
      call = NULL
    )
  }
  invisible(NULL)
}

#' Refuse an estimator that has no variance to compute inference from
#'
#' The requirement of everything that reads the variance: the standard errors
#' and everything built on them, and the influence functions, which need a bread
#' that can be inverted for the same reason the variance does.
#'
#' Two states are refused here, and the remedy for one of them does not exist
#' for the other. An estimator that was never solved has no estimates either, so
#' the caller is told to solve it, which is what `check_has_estimates()` says. A
#' fit that ran and came back without a variance has its estimates, and telling
#' that caller to fit the model names a step they already took: what happened is
#' that the bread could not be inverted, which `estimate()` warned about at the
#' time. The estimates are the property that tells the two apart, since a solve
#' records them whatever became of the sandwich, and they are what the second
#' message sends the reader to.
#'
#' @param object The estimator inference was asked for.
#' @returns Invisible `NULL` where the fit carries a variance; otherwise it
#'   throws.
#' @noRd
check_estimated <- function(object) {
  check_has_estimates(object)
  if (!is.null(object@variance)) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "This fit has no variance, so there is nothing to compute
             inference from.",
      "i" = "{.fn estimate} solved the estimating equations and recorded the
             estimates, but the bread matrix came back holding values that are
             not finite, so it could not be inverted and no sandwich variance
             was built. The fit warned about that at the time.",
      "i" = "{.fn coef} still reports the estimates. A bread that cannot be
             inverted usually means the estimating functions are not
             differentiable at the estimates, or that the finite differences
             behind it are not finite there; {.code deriv_method = \"exact\"}
             takes no difference at all."
    ),
    call = NULL
  )
}

#' @noRd
check_null_length <- function(null, n_params) {
  if (length(null) != 1L && length(null) != n_params) {
    cli::cli_abort(
      "{.arg null} must be a single value or one value per parameter
       ({n_params}), not length {length(null)}.",
      call = NULL
    )
  }
  invisible(NULL)
}

#' @noRd
check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    cli::cli_abort(
      "{.arg alpha} must be a single number between 0 and 1 (exclusive).",
      call = NULL
    )
  }
}

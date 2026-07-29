#' Broom tidiers for deli estimators
#'
#' [generics::tidy()] and [generics::glance()] methods for `MEstimator` and
#' `GMMEstimator` objects. These allow deli results to flow into tidyverse
#' pipelines. Requires the \pkg{generics} package (or \pkg{broom}) to be
#' loaded.
#'
#' @param x A fitted `MEstimator` or `GMMEstimator` object.
#' @param conf.int Logical. Include confidence intervals? Default `FALSE`.
#' @param conf.level Numeric confidence level for intervals. Default `0.95`.
#' @param ... Additional arguments (unused).
#'
#' @returns
#' - `tidy()`: A data.frame with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `s.value`. If `conf.int = TRUE`, also includes
#'   `conf.low` and `conf.high`.
#' - `glance()`: A single-row data.frame with model-level summaries.
#'
#' @seealso [deli-augment], the third broom generic, which returns the
#'   observation-level fitted values, intervals, and residuals.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' # Attaching generics or broom makes the shorter tidy(fit) and glance(fit)
#' # calls work as well
#' generics::tidy(fit, conf.int = TRUE)
#'
#' generics::glance(fit)
#'
#' @name deli-tidiers
#' @importFrom generics tidy glance
NULL

# ---- External generic declarations ------------------------------------------

generics_tidy <- new_external_generic("generics", "tidy", "x")
generics_glance <- new_external_generic("generics", "glance", "x")

# ---- tidy --------------------------------------------------------------------

method(generics_tidy, deli_estimator) <- function(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  ...
) {
  tidy_estimator(x, conf.int = conf.int, conf.level = conf.level)
}

#' @noRd
tidy_estimator <- function(x, conf.int = FALSE, conf.level = 0.95) {
  check_estimated(x)

  theta <- x@theta
  se <- sqrt(diag(x@variance))
  z <- z_scores(x)
  p <- p_values(x)
  s <- s_values(x)

  result <- data.frame(
    term = names(theta),
    estimate = unname(theta),
    std.error = unname(se),
    statistic = unname(z),
    p.value = unname(p),
    s.value = unname(s),
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    ci <- confidence_intervals(x, alpha = 1 - conf.level)
    result$conf.low <- unname(ci[, "lower"])
    result$conf.high <- unname(ci[, "upper"])
  }

  result
}

# ---- glance ------------------------------------------------------------------

method(generics_glance, deli_estimator) <- function(x, ...) {
  check_estimated(x)
  data.frame(
    nobs = x@n_obs,
    npar = x@n_params,
    estimator = S7::S7_class(x)@name,
    finite_correction = x@finite_correction %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

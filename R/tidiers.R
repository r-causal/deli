#' Broom tidiers for deli estimators
#'
#' `tidy()` and `glance()` methods for `MEstimator` and `GMMEstimator` objects.
#' These allow deli results to flow into tidyverse pipelines.
#'
#' @param x A fitted `MEstimator` or `GMMEstimator` object.
#' @param conf.int Logical. Include confidence intervals? Default `FALSE`.
#' @param conf.level Numeric confidence level for intervals. Default `0.95`.
#' @param ... Not used. `tidy()` requires them to be empty, so that a misspelled
#'   `conf.int` or `conf.level` is an error rather than a table silently
#'   returned without intervals or at the default level. `glance()` has no
#'   optional argument for a wrong name to displace and ignores them.
#'
#' @returns
#' - `tidy()`: A data.frame with columns `term`, `estimate`, `std.error`,
#'   `statistic`, `p.value`, `s.value`. If `conf.int = TRUE`, also includes
#'   `conf.low` and `conf.high`. A `p.value` that underflows to exactly zero is
#'   reported as `0` alongside an infinite `s.value`; see [s_values()].
#' - `glance()`: A single-row data.frame with model-level summaries: `nobs`,
#'   `npar`, `estimator`, `finite_correction`, and the Hansen J-statistic of an
#'   over-identified GMM fit in `j_statistic`, `j_df` and `j_p_value`. The three
#'   J columns are present on every fit and hold the typed missing value of
#'   their own type where there is no such statistic, which is every
#'   M-estimation fit and every just-identified or `subset` GMM fit; see
#'   [GMMEstimator()] for what the statistic reads and where it is left unset.
#'
#' @seealso [deli-augment], the third broom generic, which returns the
#'   observation-level fitted values, intervals, and residuals, and [reexports]
#'   for the generics themselves, which deli re-exports so that `tidy()` and
#'   `glance()` resolve with deli alone attached.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' tidy(fit, conf.int = TRUE)
#'
#' glance(fit)
#'
#' @name deli-tidiers
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
  # `sys.call(-1)` names the `tidy()` call the caller wrote; see the comment on
  # the same call in `predict()` for why the default would name the wrong frame.
  rlang::check_dots_empty(call = sys.call(-1))
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
  # Every column here is recorded by the solve, so a fit whose variance could
  # not be built has all of them. `tidy()` above reports standard errors and
  # asks for the variance.
  check_has_estimates(x)
  j <- glance_j_statistic(x)
  data.frame(
    nobs = x@n_obs,
    npar = x@n_params,
    estimator = S7::S7_class(x)@name,
    finite_correction = x@finite_correction %||% NA_character_,
    j_statistic = j$statistic,
    j_df = j$df,
    j_p_value = j$p_value,
    stringsAsFactors = FALSE
  )
}

#' The over-identification reading of a fit, for one row of `glance()`
#'
#' The J-statistic exists for an over-identified GMM fit and for nothing else:
#' an M-estimation fit has no such reading at all, and a just-identified or a
#' `subset` GMM fit leaves the property empty for the reasons [GMMEstimator()]
#' gives. Those fits still get the three columns, holding the typed missing
#' value of the column's own type, so a set of fits bound together keeps one
#' column per quantity rather than gaining and losing columns by estimator.
#'
#' The degrees of freedom and the P-value are computed here rather than stored,
#' as [`summary()`][deli-display] computes them, so the two reports cannot
#' disagree about a statistic recorded once.
#'
#' @param x A fitted estimator.
#'
#' @returns A list holding `statistic`, `df` and `p_value`, each a length-one
#'   vector of the type its column carries.
#' @noRd
glance_j_statistic <- function(x) {
  empty <- list(
    statistic = NA_real_,
    df = NA_integer_,
    p_value = NA_real_
  )
  if (!S7::S7_inherits(x, GMMEstimator) || is.null(x@j_statistic)) {
    return(empty)
  }
  # The weight matrix is indexed by the moment conditions on every path, so its
  # dimension is the count the degrees of freedom are measured from. A `subset`
  # fit never reaches here with a statistic, so the parameters below are all of
  # them.
  df <- nrow(x@weight_matrix) - length(x@theta)
  list(
    statistic = as.numeric(x@j_statistic),
    df = as.integer(df),
    p_value = stats::pchisq(x@j_statistic, df, lower.tail = FALSE)
  )
}

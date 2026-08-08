# ---- cli brace escaping -----------------------------------------------------
# cli reads a brace in a string it formats as the start of an expression to
# evaluate, and a parameter name is written by whoever wrote the estimating
# function, so it can hold anything. Doubling each brace is how cli's own parser
# is told to take one literally; cli has an internal helper that does this, but
# it is not exported.
#
# Escaping rather than validating is the only workable choice, because cli does
# not reliably refuse a name it cannot format. "{foo}" and "beta{" do raise an
# error, but "theta{1}" prints as theta1 and "E[Y^{a=1}]" as "E[Y^1: ]", each
# without a word of complaint. A name the package's own accessors return intact
# has to survive its printing too.
#
# `NULL` is passed straight back rather than escaped, because `gsub()` answers it
# with `character(0)` and `stats::setNames()` pads a short name vector with `NA`,
# so an escaped `NULL` would relabel every coefficient `NA`. Nothing hands it one:
# `print_estimator()` resolves the labels through `default_param_names()` before
# escaping them, so what arrives here always has one entry per displayed
# parameter.

#' @noRd
escape_cli_braces <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- gsub("{", "{{", x, fixed = TRUE)
  gsub("}", "}}", x, fixed = TRUE)
}

#' @noRd
print_estimator <- function(x, label, subset = NULL) {
  check_estimator_subset(subset, x@n_params)

  if (is.null(x@theta)) {
    cli::cli_text("{.cls {label}}")
    cli::cli_bullets(c(
      " " = "Parameters: {x@n_params}",
      "i" = "Call {.fn estimate} to fit."
    ))
  } else {
    cli::cli_text("{.cls {label}}")
    cli::cli_bullets(c(
      " " = "Parameters: {x@n_params}",
      " " = "Observations: {x@n_obs}"
    ))
    theta <- x@theta
    # Filled before the subset is taken, so a numbered label reports which
    # parameter of the whole fit is on show rather than its position among the
    # displayed ones, which is how the summary path labels a subset as well.
    names(theta) <- default_param_names(names(theta), length(theta))
    if (!is.null(subset)) {
      theta <- theta[subset]
    }
    theta_fmt <- format(round(theta, 4), nsmall = 4)
    cli::cli_text("Coefficients:")
    # cli_dl() interpolates its labels, so the parameter names are escaped on
    # the way in. The values are formatted numbers and carry no braces.
    cli::cli_dl(stats::setNames(theta_fmt, escape_cli_braces(names(theta))))
  }

  invisible(x)
}

#' Display methods for deli estimators
#'
#' Methods for [base::print()] and [base::summary()] so that a fitted estimator
#' shows its coefficients at the console and reports its parameter-level
#' inference as a table.
#'
#' @details
#' `print()` reports the parameter and observation counts followed by the
#' estimates, rounded to four decimal places. An estimator that has not been
#' through [estimate()] has no estimates to show, so it reports its parameter
#' count and says so.
#'
#' `summary()` collects the estimates, their standard errors, Z-scores,
#' confidence limits, P-values, and S-values into one object, which prints as a
#' table with one row per parameter. The values are those
#' [`confint()`][deli-generics], [z_scores()], [p_values()], and [s_values()]
#' return individually, so `alpha` here means what it means there: `0.05` gives
#' 95% limits.
#'
#' An over-identified `GMMEstimator` fit reports Hansen's J-statistic above the
#' table, with its degrees of freedom and its P-value, because it judges the fit
#' as a whole rather than any one parameter. See [GMMEstimator()] for what it
#' means and where its reference distribution holds. No other fit has one, so no
#' other output carries the line.
#'
#' The `S` column reads `Inf` when a P-value underflows to exactly zero, for the
#' reason [s_values()] gives. The `P` column reports that same underflow as
#' `<2e-16`: [base::format.pval()] stops printing digits below the `eps` it is
#' given, and the table gives it `2.2e-16`. [`tidy()`][deli-tidiers] returns the
#' literal `0` instead.
#'
#' `subset` restricts which parameters are displayed and nothing else. The
#' reported parameter count and every reported value are computed from the whole
#' fit, so displaying a subset of a stacked estimator is a way to read the
#' parameters of interest without the nuisance parameters, not a way to refit
#' without them. Row labels keep the names the full fit gave them.
#'
#' @param x,object A fitted `MEstimator` or `GMMEstimator` object. Named `x` for
#'   `print()` and `object` for `summary()`, because [base::print()] and
#'   [base::summary()] name their first argument that.
#' @param alpha Numeric significance level for the confidence limits reported by
#'   `summary()`, between 0 and 1. Default `0.05` for 95% limits.
#' @param subset Integer vector of parameter indices to display, or `NULL`
#'   (default) to display all of them.
#' @param ... Not used. Must be empty, so that a name neither method recognizes
#'   is an error rather than silently ignored: a misspelled `alpha` would report
#'   limits at the default width and a misspelled `subset` would display every
#'   parameter, each while the call still read as the one that was meant.
#'
#' @returns
#' - `print()`: its input, invisibly.
#' - `summary()`: an object holding the estimates, standard errors, Z-scores,
#'   confidence limits, P-values, and S-values, which prints as a table.
#'
#' @seealso [deli-generics] for the accessors that return these quantities as
#'   plain vectors and matrices, and [deli-tidiers] for the same results as a
#'   data frame.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' fit
#'
#' # `subset` is an argument of the method rather than of the fit, so showing it
#' # here means calling the generic by name.
#' print(fit, subset = 2:3)
#'
#' summary(fit)
#'
#' # `summary()` takes `subset` as well, alongside `alpha` for the width of the
#' # reported limits.
#' summary(fit, subset = 2:3, alpha = 0.1)
#'
#' @name deli-display
NULL

# The S3 class of an S7 object is package-qualified, so the methods below are
# bound to backticked names and registered through NAMESPACE S3method
# directives; see the comment above the accessors in R/generics-s3.R.

#' @rdname deli-display
#' @export
`print.deli::deli_estimator` <- function(x, ..., subset = NULL) {
  # `sys.call(-1)` names the `print()` call the caller wrote; see the comment on
  # the same call in `predict()` for why the default would name the wrong frame.
  # `subset` follows the dots, so it matches only exactly and every misspelling
  # of it reaches this guard.
  rlang::check_dots_empty(call = sys.call(-1))
  print_estimator(x, S7::S7_class(x)@name, subset = subset)
}

#' Summary results from an estimator
#'
#' S7 class holding inference results: point estimates, standard errors,
#' confidence intervals, Z-scores, P-values, and S-values.
#'
#' @noRd
EstimatorSummary <- new_class(
  "EstimatorSummary",
  properties = list(
    theta = class_numeric,
    se = class_numeric,
    ci = class_any,
    z = class_numeric,
    p = class_numeric,
    s = class_numeric,
    n_obs = class_integer,
    n_params = class_integer,
    alpha = class_numeric,
    estimator = class_character,
    j_statistic = new_property(
      class = NULL | class_double,
      default = NULL
    ),
    j_df = new_property(
      class = NULL | class_integer,
      default = NULL
    )
  )
)

#' @noRd
summarize_estimator <- function(
  object,
  estimator_label,
  alpha = 0.05,
  subset = NULL
) {
  check_estimated(object)
  check_alpha(alpha)
  check_estimator_subset(subset, object@n_params)

  theta <- object@theta
  se <- sqrt(diag(object@variance))
  ci <- confidence_intervals(object, alpha = alpha)
  z <- z_scores(object)
  p <- p_values(object)
  s <- s_values(object)

  # Only a GMMEstimator carries a J-statistic, and only an over-identified fit
  # has one set; see GMMEstimator() for what it reports. Its degrees of freedom
  # are the moment conditions the weight matrix is indexed by, beyond the
  # parameters. A `subset` fit never reaches here with one, so the count below is
  # the number of parameters the fit estimated.
  j_statistic <- NULL
  j_df <- NULL
  if (S7::S7_inherits(object, GMMEstimator) && !is.null(object@j_statistic)) {
    j_statistic <- object@j_statistic
    j_df <- nrow(object@weight_matrix) - length(object@theta)
  }

  # A subset restricts the displayed rows to the given parameter indices while
  # leaving the reported values and n_params untouched. The row labels track the
  # original parameter index because the named theta vector carries its names
  # through the subset.
  if (!is.null(subset)) {
    names(theta) <- default_param_names(names(theta), length(theta))
    theta <- theta[subset]
    se <- se[subset]
    ci <- ci[subset, , drop = FALSE]
    z <- z[subset]
    p <- p[subset]
    s <- s[subset]
  }

  EstimatorSummary(
    theta = theta,
    se = se,
    ci = ci,
    z = z,
    p = p,
    s = s,
    n_obs = object@n_obs,
    n_params = object@n_params,
    alpha = alpha,
    estimator = estimator_label,
    j_statistic = j_statistic,
    j_df = j_df
  )
}

#' @rdname deli-display
#' @export
`summary.deli::deli_estimator` <- function(
  object,
  alpha = 0.05,
  subset = NULL,
  ...
) {
  rlang::check_dots_empty(call = sys.call(-1))
  summarize_estimator(
    object,
    S7::S7_class(object)@name,
    alpha = alpha,
    subset = subset
  )
}

#' @export
`print.deli::EstimatorSummary` <- function(x, ...) {
  cli::cli_rule("{x@estimator} Results")

  cli::cli_text("Observations: {x@n_obs}")
  cli::cli_text("Parameters:   {x@n_params}")
  # An over-identified GMM fit reports the J-statistic of its moment conditions
  # alongside the counts, since it describes the fit as a whole rather than any
  # one parameter. The P-value is formatted as the table's own column is, so the
  # two agree on how small a probability is reported.
  #
  # Both halves are required, rather than the statistic alone, because the
  # degrees of freedom are what the statistic is read against: summarize_estimator()
  # fills the two together, but the class can be constructed directly, and a
  # statistic with no degrees of freedom would reach pchisq() as a missing
  # argument rather than being left unreported.
  if (!is.null(x@j_statistic) && !is.null(x@j_df)) {
    j_fmt <- format(round(x@j_statistic, 4), nsmall = 4)
    j_p_fmt <- format.pval(
      stats::pchisq(x@j_statistic, x@j_df, lower.tail = FALSE),
      digits = 3,
      eps = 2.2e-16
    )
    cli::cli_text("J-statistic:  {j_fmt} on {x@j_df} df (P = {j_p_fmt})")
  }
  cli::cli_text("")

  # Build the results table
  pct <- (1 - x@alpha) * 100
  lcl_lab <- sprintf("%.0f%% LCL", pct)
  ucl_lab <- sprintf("%.0f%% UCL", pct)

  # Determine column widths. The rows below go out through cli_verbatim(),
  # which prints its argument as given, so the names need no escaping here;
  # print_estimator() escapes because cli_dl() interpolates its labels.
  param_names <- default_param_names(names(x@theta), length(x@theta))
  name_width <- max(nchar(param_names), 5)

  # Header
  header <- sprintf(
    "%-*s %10s %10s %10s %10s %10s %10s %10s",
    name_width,
    "",
    "Estimate",
    "Std.Err",
    "Z-score",
    lcl_lab,
    ucl_lab,
    "P-value",
    "S-value"
  )
  cli::cli_verbatim(header)

  # Rows
  for (i in seq_along(x@theta)) {
    p_fmt <- format.pval(x@p[i], digits = 3, eps = 2.2e-16)
    row <- sprintf(
      "%-*s %10.4f %10.4f %10.4f %10.4f %10.4f %10s %10.4f",
      name_width,
      param_names[i],
      x@theta[i],
      x@se[i],
      x@z[i],
      x@ci[i, "lower"],
      x@ci[i, "upper"],
      p_fmt,
      x@s[i]
    )
    cli::cli_verbatim(row)
  }

  invisible(x)
}

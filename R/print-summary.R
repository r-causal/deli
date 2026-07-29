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
# `NULL` is passed straight back rather than escaped. `gsub()` answers it with
# `character(0)`, and `stats::setNames()` pads a short name vector with `NA`, so
# escaping an unnamed theta would relabel every coefficient `NA` instead of
# leaving cli to refuse an unnamed list. A fitted estimator always names its
# parameters; an unnamed one is reachable only by assigning `@theta` directly,
# and that should keep failing loudly.

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

#' @noRd
method(print, deli_estimator) <- function(x, ..., subset = NULL) {
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
    estimator = class_character
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
    estimator = estimator_label
  )
}

#' @noRd
method(summary, deli_estimator) <- function(
  object,
  alpha = 0.05,
  subset = NULL,
  ...
) {
  summarize_estimator(
    object,
    S7::S7_class(object)@name,
    alpha = alpha,
    subset = subset
  )
}

#' @noRd
method(print, EstimatorSummary) <- function(x, ...) {
  cli::cli_rule("{x@estimator} Results")

  cli::cli_text("Observations: {x@n_obs}")
  cli::cli_text("Parameters:   {x@n_params}")
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

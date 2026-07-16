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
    cli::cli_dl(stats::setNames(theta_fmt, names(theta)))
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
summarize_estimator <- function(object, estimator_label, alpha = 0.05,
                                subset = NULL) {
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
    names(theta) <- names(theta) %||% paste0("theta_", seq_along(theta))
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
method(summary, deli_estimator) <- function(object, alpha = 0.05,
                                            subset = NULL, ...) {
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

  # Determine column widths
  param_names <- names(x@theta)
  if (is.null(param_names)) {
    param_names <- paste0("theta_", seq_along(x@theta))
  }
  name_width <- max(nchar(param_names), 5)

  # Header
  header <- sprintf(
    "%-*s %10s %10s %10s %10s %10s %10s %10s",
    name_width, "", "Estimate", "Std.Err", "Z-score",
    lcl_lab, ucl_lab, "P-value", "S-value"
  )
  cli::cli_verbatim(header)

  # Rows
  for (i in seq_along(x@theta)) {
    p_fmt <- format.pval(x@p[i], digits = 3, eps = 2.2e-16)
    row <- sprintf(
      "%-*s %10.4f %10.4f %10.4f %10.4f %10.4f %10s %10.4f",
      name_width, param_names[i],
      x@theta[i], x@se[i], x@z[i],
      x@ci[i, "lower"], x@ci[i, "upper"],
      p_fmt, x@s[i]
    )
    cli::cli_verbatim(row)
  }

  invisible(x)
}

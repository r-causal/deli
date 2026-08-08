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

# The S3 class of an S7 object is package-qualified, so the methods below are
# bound to backticked names and registered through NAMESPACE S3method
# directives; see the comment above the accessors in R/generics-s3.R.

# ---- tidy --------------------------------------------------------------------

#' @rdname deli-tidiers
#' @export
`tidy.deli::deli_estimator` <- function(
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

#' @rdname deli-tidiers
#' @export
`glance.deli::deli_estimator` <- function(x, ...) {
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

# ---- What augment() decides that predict() does not --------------------------
# `augment()` computes nothing of its own. Every column it adds comes from
# `predict()` or `residuals()`, and the whole of its design is three choices
# about how to present them.
#
# The first is the scale of `.fitted`. It is the link scale by default, with a
# `type.predict` argument to ask for the response scale, which is what
# `broom::augment()` does for a `glm` and what this package's own `predict()`
# already defaults to. Following both leaves one rule to remember rather than
# two, and makes `.fitted` equal to `predict()` on the same fit under the
# defaults of each.
#
# The argument is spelled `type.predict` rather than `type` for the same reason
# `predict()` spells its arguments `se.fit` and `level`: it is the name the
# ecosystem uses at this entry point. broom's augment methods take
# `type.predict` because an augmented frame carries fitted values and residuals
# at once, so a bare `type` would not say which it governed. That is true here
# as well.
#
# The second is that `.resid` does not follow `type.predict`. It is always the
# response residual `residuals()` returns, the response minus its conditional
# mean, because a residual against the linear predictor would subtract a
# log-odds from a 0/1 response. broom pairs a link-scale `.fitted` with a
# deviance residual for a `glm` and so makes the same split; this package cannot
# offer a deviance residual, since an M-estimator need not come from a
# likelihood, so the response residual is what stands in every column.
#
# The third is that a `newdata` frame comes back without `.resid` at all, rather
# than with a column of `NA`. New covariate values carry no response in general,
# and a column that is present for some data and absent for other data is worse
# to program against than one that is absent whenever it cannot be computed.

#' Augment data with predictions from a fitted deli estimator
#'
#' An `augment()` method for `MEstimator` and `GMMEstimator` objects fitted
#' through the formula interface. It returns the model frame the fit was built
#' from, or `newdata` when supplied, with the fitted values, their standard
#' errors, a Wald confidence interval, and the residuals as columns beside it.
#'
#' @details
#' The columns added are `.fitted`, `.se.fit`, `.lower`, `.upper`, and, when
#' `newdata` is not supplied, `.resid`. They are exactly what
#' [`predict()`][deli-predict] and [`residuals()`][deli-generics] return for the
#' same fit, so `.fitted` is `predict()`, `.lower` and `.upper` are
#' `predict(interval = "confidence")`, and `.resid` is `residuals()`.
#'
#' `.fitted` is on the link scale by default, matching both
#' [`predict()`][deli-predict] and `broom::augment()` on a `glm`, and
#' `type.predict = "response"` puts it and its interval on the scale of the
#' response. `.resid` is the response residual either way, since a residual
#' measured against a linear predictor would be a response minus a quantity the
#' response is not measured in.
#'
#' Rows the fit dropped for missing data are not reported, so the result has one
#' row per [`nobs()`][deli-generics] and its row names are those of the retained
#' rows. A `newdata` row with a missing value is kept, with `NA` in the added
#' columns, so that the result lines up with the rows handed in.
#'
#' `augment()` covers the estimating equations whose linear predictor
#' [`predict()`][deli-predict] forms, and refuses the same fits with the same
#' reasons; see [deli-predict]. It has no counterpart to that method's `times`
#' argument: a survival measure is one value per row of the data and time rather
#' than one per row, so the predictions do not go beside the data as columns.
#'
#' @param x A fitted `MEstimator` or `GMMEstimator` object made with the formula
#'   interface (after calling [estimate()]).
#' @param newdata A data frame of covariate values to predict at, or `NULL`
#'   (default) to report the model frame the fit was built from. That frame
#'   holds the variables the formula named, so a transformed term appears as the
#'   column the transformation produced and a column of the fitting data the
#'   formula did not name is not there.
#' @param type.predict Character string. `"link"` (default) puts `.fitted` and
#'   its interval on the scale of the linear predictor; `"response"` puts them
#'   on the scale of the response.
#' @param conf.level Numeric confidence level for `.lower` and `.upper`.
#'   Default `0.95`.
#' @param ... Not used. Must be empty, so that a name that is not one of the
#'   documented arguments is an error rather than silently ignored. A
#'   misspelled `newdata` would otherwise augment the fitted rows while the
#'   caller believed they had asked for rows of their own.
#'
#' @returns A data frame: the model frame, or `newdata`, followed by the columns
#'   `.fitted`, `.se.fit`, `.lower`, `.upper`, and `.resid`. The last is absent
#'   when `newdata` is supplied.
#'
#' @seealso [deli-predict] for the predictions themselves and the equations they
#'   are available for, [deli-tidiers] for the parameter-level and model-level
#'   summaries, and [reexports] for the generic itself, which deli re-exports so
#'   that `augment()` resolves with deli alone attached.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' head(augment(fit))
#'
#' # New covariate patterns come back with no residual column, since they carry
#' # no response to residualize against.
#' augment(fit, newdata = data.frame(wt = c(2, 3, 4), hp = 110))
#'
#' @name deli-augment
NULL

# The S3 class of an S7 object is package-qualified, so the method below is
# bound to a backticked name and registered through a NAMESPACE S3method
# directive; see the comment above the accessors in R/generics-s3.R.

# ---- augment -----------------------------------------------------------------

#' @rdname deli-augment
#' @export
`augment.deli::deli_estimator` <- function(
  x,
  newdata = NULL,
  type.predict = c("link", "response"),
  conf.level = 0.95,
  ...
) {
  # `sys.call(-1)` names the `augment()` call the caller wrote; see the comment
  # on the same call in `predict()` for why the default would name the wrong
  # frame.
  rlang::check_dots_empty(call = sys.call(-1))
  check_estimated(x)
  type.predict <- match.arg(type.predict)
  check_level(conf.level, "conf.level")
  check_augment_newdata(newdata)
  resolve_predict_link(x, "augment")

  prediction <- predict(
    x,
    newdata = newdata,
    type = type.predict,
    se.fit = TRUE,
    interval = "confidence",
    level = conf.level
  )

  augmented <- augment_frame(x@model_spec, newdata)
  augmented$.fitted <- unname(prediction$fit[, "fit"])
  augmented$.se.fit <- unname(prediction$se.fit)
  augmented$.lower <- unname(prediction$fit[, "lwr"])
  augmented$.upper <- unname(prediction$fit[, "upr"])
  if (is.null(newdata)) {
    augmented$.resid <- unname(residuals(x))
  }
  augmented
}

# ---- The data the columns are added to ---------------------------------------

#' The frame `augment()` reports against
#'
#' `newdata` when there is one, and otherwise the model frame the fit was built
#' from, stripped of the `terms` and `na.action` attributes that make it a model
#' frame. What comes back is a table of values rather than a frame anything else
#' is built from, and those attributes would travel with it into whatever the
#' caller does next.
#'
#' @param spec The fit's `model_spec`.
#' @param newdata A data frame, or `NULL`.
#' @returns A data frame.
#' @noRd
augment_frame <- function(spec, newdata) {
  if (!is.null(newdata)) {
    return(newdata)
  }
  mf <- spec$model_frame
  attr(mf, "terms") <- NULL
  attr(mf, "na.action") <- NULL
  mf
}

# ---- Validation --------------------------------------------------------------

#' @noRd
check_augment_newdata <- function(newdata) {
  if (!is.null(newdata) && !is.data.frame(newdata)) {
    cli::cli_abort(
      c(
        "{.arg newdata} must be a data frame or {.val NULL}, not
         {.obj_type_friendly {newdata}}.",
        "i" = "{.fn augment} returns it with the fitted columns beside its own
               columns, so it has to be something columns can be added to."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

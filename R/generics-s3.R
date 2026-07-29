# ---- What the five model accessors report on ---------------------------------
# `fitted()`, `residuals()`, `model.frame()`, `formula()`, and `terms()` all
# report on the model a formula fit was specified as, which only
# `formula_model_spec()` records. A fit built from a `stacked_equations` closure
# has none, so each of them is an error there, naming itself and saying why.
#
# `fitted()` is the conditional mean of the response, so it is
# `predict(type = "response")` and not the linear predictor. That is the
# convention of `fitted.glm()`, whose fitted values are on the response scale
# while `predict.glm()` answers on the link scale by default, and it is the only
# reading under which `residuals()` is the response minus the fitted value.
#
# `residuals()` is that difference and nothing else. Under a non-identity link
# it is what a GLM calls a response residual, `residuals(type = "response")`:
# the response minus E[Y | X], both on the scale the response was measured in.
# The alternatives a GLM offers are not available here. A deviance residual
# needs a likelihood, and an M-estimator need not come from one, so there is
# none to take for a robust or penalized fit; a Pearson residual needs a
# variance function, which the estimating equation does not have to state
# either. The response residual is the one definable for every equation
# `predict()` supports, and it is defined the same way for all of them, so
# nothing about the link changes what the subtraction means. What the link does
# change is the spread: a response residual is heteroscedastic by construction
# under a non-identity link, and is not the quantity to compare across
# observations without standardizing it.
#
# `model.frame()`, `formula()`, and `terms()` are the specification rather than
# the solve, and do not call `check_estimated()`. Nothing they return depends on
# the estimates, and reporting "cannot compute inference before calling
# estimate()" from `formula()` would name a step the caller did not ask for. A
# fit that carries a specification but no estimates cannot arise through the
# formula interface anyway, which estimates in the same call that records it.
#
# The three of them take their first argument under the name the base generic
# dispatches on: `x` for `formula()` and `terms()`, and `formula` for
# `model.frame()`. That last one reads as though it takes a formula and does
# not; the name is the generic's and dispatch follows it.
#
# ---- Why model.matrix() has a method of its own ------------------------------
# `model.matrix.default()` would answer for a fit with no method written at all:
# it reads the terms and the model frame off the object through the two generics
# above and builds from them. What it builds the factor columns with is
# `getOption("contrasts")` rather than the coding the fit recorded, so a fit made
# under sum contrasts and read back in a session under the default treatment
# contrasts comes back as a design of columns the fit never had, and nothing in
# the result says the coding changed. A design matrix is the fit's rather than the
# session's, which is the rule `predict()` already follows; both go through
# `rebuild_design()`, so the design a fit reports is the design its predictions
# are formed on.
#
# Without `data` the recorded design is reported as recorded, which is the matrix
# a rebuild from the recorded frame produces anyway. With `data`, the response is
# deleted from the terms first, so covariate values that carry none are enough.
# That is where a design parts company with a model frame, which holds every
# variable the formula names and says so when the response is absent. There is no
# `contrasts.arg`: the argument that would take one is the fit, which has already
# answered the question.
#
# ---- Why model.frame() takes the rest of its formals -------------------------
# `model.frame()` honors `data` for the reason every base method does. A method
# that swallowed it into `...` would answer with the frame the fit was solved on,
# at the fitted number of rows, and give no sign that the caller's rows had been
# discarded.
#
# With `data`, the frame is rebuilt through the recorded terms and `xlevels`,
# the same triple `predict()` rebuilds a design from and for the same reasons:
# the terms carry the `predvars` of a term such as `poly(x, 2)`, so its basis is
# the fitted one rather than one refitted to the rows handed in, and the fitted
# factor levels give a factor covering only some of them the fitted set of
# columns. `subset` and `na.action` are passed on to `stats::model.frame()` and
# mean there what they mean everywhere else, including the `na.omit` default,
# which is where this differs from `predict()`: a prediction lines up with the
# rows of `newdata`, and a model frame is a frame built under an `na.action`.
#
# The response has to be in `data`. A model frame holds every variable the
# formula names, the response included, which is why `model.frame()` on an `lm`
# fit requires it as well. Covariate values that carry no response are what
# `predict()` and `augment()` take a `newdata` for, and the error says so.
#
# Without `data` there is nothing to build, so the recorded frame comes back as
# recorded and an argument that describes how to build one is an error rather
# than a silent no-op. `xlev` is what makes that rule delicate: `NULL` is its own
# default rather than something a caller asks for, so a `NULL` `xlev` beside a
# `NULL` `data` has to stay legal and only a supplied one is refused.

#' Standard S3 generics for deli estimators
#'
#' Methods for base R generics [stats::coef()], [stats::vcov()],
#' [stats::confint()], [stats::nobs()], [stats::fitted()],
#' [stats::residuals()], [stats::model.frame()], [stats::model.matrix()],
#' [stats::formula()], and [stats::terms()] so that deli estimator objects
#' interoperate with the broader R modeling ecosystem.
#'
#' @details
#' The last six report on the model a fit was specified as, which only the
#' formula interface records, so each of them is an error for a fit built from a
#' `stacked_equations` function. `fitted()` and `residuals()` are also an error
#' for a formula fit of an estimating equation [`predict()`][deli-predict] does
#' not support, since a fitted value is a prediction; see [deli-predict] for the
#' equations it covers.
#'
#' `fitted()` is the conditional mean of the response, so it is on the response
#' scale rather than the link scale, matching [stats::fitted()] on a `glm`
#' object. `residuals()` is the response minus that mean, the residual a GLM
#' calls a response residual, and the only `type` it takes. Deviance and Pearson
#' residuals are not offered: the first needs a likelihood and the second a
#' variance function, and an M-estimator need not state either, while the
#' response residual is defined the same way for every equation `predict()`
#' supports. Under a non-identity link it is heteroscedastic by construction.
#'
#' `model.frame()` returns the frame the fit was solved on, or builds the frame
#' of `data` when one is supplied, through the terms and factor levels the fit
#' recorded. A data-dependent term such as `poly(x, 2)` is therefore evaluated
#' at its fitted basis and a factor gets its fitted levels, as they are under
#' [`predict()`][deli-predict]. `data` has to carry the response, since a model
#' frame holds every variable the formula names; covariate values that carry no
#' response are what `predict()` and [deli-augment] take a `newdata` for.
#'
#' `model.matrix()` returns the design the fit was solved on, or the design of
#' `data` when one is supplied, coded with the contrasts and factor levels the fit
#' recorded rather than with whatever `getOption("contrasts")` says when the call
#' is made. A design matrix is a property of the fit, so a fit made under one
#' setting of that option answers with the coding it was solved on under any
#' other, which is what [`predict()`][deli-predict] does as well. A design names
#' no response, so `data` needs only the predictors.
#'
#' @param object A fitted `MEstimator` or `GMMEstimator` object (after calling
#'   [estimate()]).
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level The confidence level required. Default `0.95`.
#' @param type Character string naming the residual `residuals()` returns. Only
#'   `"response"` is available, since it is the one residual defined for every
#'   estimating equation `predict()` supports. Any other value is an error
#'   rather than a response residual under another name.
#' @param x A fitted `MEstimator` or `GMMEstimator` object. Named `x` because
#'   [stats::formula()] and [stats::terms()] name their first argument that.
#' @param formula A fitted `MEstimator` or `GMMEstimator` object. Named
#'   `formula` because [stats::model.frame()] names its first argument that;
#'   the method takes a fit, not a formula.
#' @param data A data frame to build the model frame or the design matrix of, or
#'   `NULL` (default) to report the one the fit was solved on. `model.frame()`
#'   needs it to carry the response as well as the predictors; `model.matrix()`
#'   names no response, so the predictors alone are enough there.
#' @param subset,na.action,drop.unused.levels,xlev Passed to
#'   [stats::model.frame()], and meaning there what they mean for any other
#'   model frame. `xlev` defaults to the factor levels the fit recorded. All
#'   four describe how to build a frame from `data`, so supplying one without
#'   `data` is an error rather than a silent no-op.
#' @param ... Additional arguments. `model.frame()` requires them to be empty,
#'   so that a name it does not recognize is an error rather than silently
#'   ignored; the rest of the methods ignore them.
#'
#' @returns
#' - `coef()`: Named numeric vector of parameter estimates.
#' - `vcov()`: Named variance-covariance matrix.
#' - `confint()`: Matrix with columns `"lower"` and `"upper"`.
#' - `nobs()`: Integer number of observations.
#' - `fitted()`: Named numeric vector of fitted values on the response scale,
#'   one per observation the fit was solved on.
#' - `residuals()`: Named numeric vector of response residuals.
#' - `model.frame()`: The model frame the fit was built from, with the rows
#'   dropped for missing data already removed, or the model frame of `data`
#'   when one is supplied.
#' - `model.matrix()`: The design matrix the fit was solved on, or the design
#'   matrix of `data` when one is supplied, coded with the contrasts and factor
#'   levels the fit recorded.
#' - `formula()`: The model formula.
#' - `terms()`: The `terms` object of the model frame, carrying the response
#'   index, any offset, and the `predvars` of a data-dependent term.
#'
#' @seealso [deli-predict] for predictions at new covariate values, and
#'   [deli-augment] for the fitted values, intervals, and residuals as columns
#'   beside the data.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' coef(fit)
#'
#' vcov(fit)
#'
#' confint(fit)
#'
#' nobs(fit)
#'
#' head(fitted(fit))
#'
#' head(residuals(fit))
#'
#' formula(fit)
#'
#' @name deli-generics
#' @importFrom stats coef vcov confint nobs fitted residuals model.frame
#' @importFrom stats model.matrix formula terms
NULL

# ---- External generic declarations ------------------------------------------

stats_coef <- new_external_generic("stats", "coef", "object")
stats_vcov <- new_external_generic("stats", "vcov", "object")
stats_confint <- new_external_generic("stats", "confint", "object")
stats_nobs <- new_external_generic("stats", "nobs", "object")
stats_fitted <- new_external_generic("stats", "fitted", "object")
stats_residuals <- new_external_generic("stats", "residuals", "object")
stats_model_frame <- new_external_generic("stats", "model.frame", "formula")
stats_model_matrix <- new_external_generic("stats", "model.matrix", "object")
stats_formula <- new_external_generic("stats", "formula", "x")
stats_terms <- new_external_generic("stats", "terms", "x")

# ---- coef --------------------------------------------------------------------

method(stats_coef, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@theta
}

# ---- vcov --------------------------------------------------------------------

method(stats_vcov, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@variance
}

# ---- confint -----------------------------------------------------------------

method(stats_confint, deli_estimator) <- function(
  object,
  parm,
  level = 0.95,
  ...
) {
  check_estimated(object)
  ci <- confidence_intervals(object, alpha = 1 - level)
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}

# ---- nobs --------------------------------------------------------------------

method(stats_nobs, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  object@n_obs
}

# ---- fitted ------------------------------------------------------------------

method(stats_fitted, deli_estimator) <- function(object, ...) {
  check_estimated(object)
  # Resolved before predicting so that a fit no prediction can be made from is
  # reported against the function the user called.
  resolve_predict_link(object, "fitted")
  predict(object, type = "response")
}

# ---- residuals ---------------------------------------------------------------

method(stats_residuals, deli_estimator) <- function(
  object,
  type = "response",
  ...
) {
  check_estimated(object)
  check_residual_type(type)
  resolve_predict_link(object, "residuals")
  mu <- predict(object, type = "response")

  # The recorded response is the vector the estimating equation was handed, so
  # it is already coded the way the fit solved for it and already filtered to
  # the rows the design has. The names are taken from the fitted values rather
  # than from it, because the coercion of a factor response to a 0/1 indicator
  # drops the ones the model frame gave it.
  residual <- as.numeric(object@model_spec$y) - unname(mu)
  names(residual) <- names(mu)
  residual
}

# ---- model.frame -------------------------------------------------------------

method(stats_model_frame, deli_estimator) <- function(
  formula,
  data = NULL,
  subset = NULL,
  na.action,
  drop.unused.levels = FALSE,
  xlev = NULL,
  ...
) {
  # `sys.call(-1)` names the `model.frame()` call the caller wrote; see the
  # comment on the same call in `predict()` for why the default would name the
  # wrong frame.
  rlang::check_dots_empty(call = sys.call(-1))
  spec <- model_spec_or_abort(formula, "model.frame")

  if (is.null(data)) {
    check_model_frame_needs_data(
      subset = subset,
      na_action = !missing(na.action),
      drop.unused.levels = drop.unused.levels,
      xlev = xlev
    )
    return(spec$model_frame)
  }

  model_terms <- spec$terms
  model_xlev <- xlev %||% spec$xlevels
  check_model_frame_response(model_terms, data)

  # The call is assembled rather than written out so that `subset` reaches
  # `stats::model.frame()` as the expression the caller wrote. It is evaluated
  # there in `data` first and in the terms' environment second, which is what
  # makes `subset = g == "a"` find the column `g`; forwarding the value this
  # frame would see instead would have evaluated it here, where there is no
  # `g`. `na.action` is passed only when the caller supplied one, so that the
  # `getOption("na.action")` default stays in the hands of the function that
  # documents it.
  args <- list(
    formula = quote(model_terms),
    data = quote(data),
    drop.unused.levels = quote(drop.unused.levels),
    xlev = quote(model_xlev)
  )
  subset_expr <- substitute(subset)
  if (!is.null(subset_expr)) {
    args$subset <- subset_expr
  }
  if (!missing(na.action)) {
    args$na.action <- quote(na.action)
  }
  eval(as.call(c(list(quote(stats::model.frame)), args)))
}

# ---- model.matrix ------------------------------------------------------------

method(stats_model_matrix, deli_estimator) <- function(
  object,
  data = NULL,
  ...
) {
  spec <- model_spec_or_abort(object, "model.matrix")

  if (is.null(data)) {
    return(spec$X)
  }

  rebuild_design(spec, data)$X
}

# ---- formula -----------------------------------------------------------------

method(stats_formula, deli_estimator) <- function(x, ...) {
  # Read off the terms rather than recorded a second time, which is how
  # `formula.lm()` answers as well, so the two records cannot disagree.
  stats::formula(model_spec_or_abort(x, "formula")$terms)
}

# ---- terms -------------------------------------------------------------------

method(stats_terms, deli_estimator) <- function(x, ...) {
  model_spec_or_abort(x, "terms")$terms
}

# ---- Validation --------------------------------------------------------------

#' Validate the residual type asked for
#'
#' @param type The `type` supplied to `residuals()`.
#' @returns Invisible `NULL`. Raises an error unless `type` is `"response"`.
#' @noRd
check_residual_type <- function(type) {
  if (identical(type, "response")) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "{.arg type} must be {.val response}.",
      "i" = "{.fn residuals} returns the response residual, the response minus
             its conditional mean, which is the one residual defined for every
             estimating equation {.fn predict} supports.",
      "i" = "The {.val deviance} and {.val pearson} residuals of a GLM are not
             available: a deviance residual needs a likelihood and a Pearson
             residual needs a variance function, and an M-estimator need not
             state either."
    ),
    call = NULL
  )
}

#' Refuse an argument that describes how to build a model frame from data
#'
#' The recorded model frame is reported as it was recorded, so an argument that
#' would have shaped a frame built from `data` has nothing to act on and is
#' refused rather than dropped.
#'
#' @param subset,drop.unused.levels,xlev The values supplied.
#' @param na_action Whether `na.action` was supplied.
#' @returns Invisible `NULL`. Raises an error if any of them was supplied.
#' @noRd
check_model_frame_needs_data <- function(
  subset,
  na_action,
  drop.unused.levels,
  xlev
) {
  supplied <- c(
    if (!is.null(subset)) "subset",
    if (na_action) "na.action",
    if (!identical(drop.unused.levels, FALSE)) "drop.unused.levels",
    if (!is.null(xlev)) "xlev"
  )
  if (length(supplied) == 0L) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "{.arg {supplied}} {?needs/need} {.arg data}.",
      "i" = "With no {.arg data}, {.fn model.frame} reports the frame the fit
             was solved on exactly as it was recorded, so there is nothing for
             {.arg {supplied}} to shape.",
      "i" = "Pass the data to build a frame from, or subset the frame that comes
             back."
    ),
    call = NULL
  )
}

#' Check that data carries the response the terms name
#'
#' @param model_terms The fit's recorded terms.
#' @param data The data a model frame is being built from.
#' @returns Invisible `NULL`. Raises an error if a response variable is absent.
#' @noRd
check_model_frame_response <- function(model_terms, data) {
  if (!is.data.frame(data) && !is.list(data)) {
    return(invisible(NULL))
  }
  response <- attr(model_terms, "response")
  if (is.null(response) || response == 0L) {
    return(invisible(NULL))
  }
  wanted <- all.vars(attr(model_terms, "variables")[[1L + response]])
  absent <- setdiff(wanted, names(data))
  if (length(absent) == 0L) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "{.arg data} does not carry the response variable{?s} {.var {absent}}.",
      "i" = "A model frame holds every variable the formula names, the response
             included, so {.fn model.frame} needs it here just as it does on an
             {.cls lm} fit.",
      "i" = "For covariate values that carry no response, use {.fn predict} or
             {.fn augment} with {.arg newdata}."
    ),
    call = NULL
  )
}

# ---- Conditions --------------------------------------------------------------

#' The model specification of a fit, or an error naming what has none
#'
#' @param object A fitted estimator.
#' @param fn The name of the function the user called, which the error names.
#' @returns The fit's `model_spec` list.
#' @noRd
model_spec_or_abort <- function(object, fn) {
  spec <- object@model_spec
  if (is.null(spec)) {
    abort_accessor_no_model_spec(fn)
  }
  spec
}

#' @noRd
abort_accessor_no_model_spec <- function(fn) {
  cli::cli_abort(
    c(
      "{.fn {fn}} needs a fit made through the formula interface.",
      "i" = "This fit was built from a {.arg stacked_equations} function, which
             records no formula, terms, or model frame to report.",
      "i" = "Pass a formula and {.arg data} to {.fn m_estimate} or
             {.fn gmm_estimate} for a fit that describes itself."
    ),
    call = NULL
  )
}

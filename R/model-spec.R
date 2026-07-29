# ---- What a model specification is for ---------------------------------------
# A formula fit knows things about itself that a fit built from a bare
# `stacked_equations` closure never can: which formula was written, how its
# factors were coded, which estimating equation was driven, and what the design
# and response looked like. None of that survives into the fitted object on its
# own, because the closure the formula interface builds is opaque to everything
# downstream of it. `formula_model_spec()` writes it down.
#
# The split between model specification and per-observation data is the
# load-bearing part. An argument such as `distribution` or `link` describes the
# model and applies unchanged to any data the model is later evaluated on. An
# argument such as `weights`, `offset`, or `event` is one value per observation
# of the data the model was fitted to, and has no counterpart in data the model
# has not seen. Only the first kind belongs in a reusable specification, so
# `spec_ee_args()` keeps those and drops the rest. The per-observation values are
# not lost: `X`, `y`, and `offset` are recorded as their own fields, and the rest
# stay reachable through the closure's environment.

#' The arguments of the built-in estimating equations that vary by observation
#'
#' Names the arguments that carry one value per row of the data rather than a
#' description of the model. Every such argument of every `ee_*()` function is
#' listed, whether or not the formula interface can reach it today, so that the
#' answer does not change when a new equation becomes reachable.
#'
#' The table is written out rather than derived, for the same reason
#' `appended_param_name()` is: the property being asked about is not visible in
#' what an estimating equation returns. Two arguments of the same length can be a
#' vector of censoring indicators and a length-2 vector of truncation bounds, and
#' only the equation's documentation says which is which. Deriving the answer
#' from the length of a value would also make the classification depend on the
#' number of observations, so one model specification would classify two ways at
#' two sample sizes.
#'
#' An argument not on the list is treated as model specification. That direction
#' is deliberate. A specification argument wrongly dropped would leave a later
#' caller silently applying a default link or distribution to new data, while a
#' per-observation argument wrongly kept is passed on at the fitted length and
#' the estimating equation rejects it. The second failure is the one a user can
#' see.
#'
#' @returns A character vector of argument names.
#' @noRd
per_observation_ee_args <- function() {
  c(
    # Shared across the regression, GLM, and survival families.
    "weights",
    "offset",
    "event",
    "time",
    "y",
    # Causal designs: treatment, instrument, and the nuisance and
    # counterfactual design matrices.
    "A",
    "Z",
    "W",
    "V",
    "X1",
    "X0",
    # Measurement error: mismeasured values, gold-standard values, and the
    # sample indicator.
    "y_star",
    "a",
    "a_star",
    "r",
    # Missingness and dose-response data.
    "delta",
    "q_eval",
    "dose",
    "response"
  )
}

#' The model-specification subset of the arguments forwarded to an equation
#'
#' Drops the per-observation arguments named by `per_observation_ee_args()` and
#' keeps the rest, which describe the model rather than the sample it was fitted
#' to.
#'
#' @param ee_args The evaluated `...` arguments forwarded to the estimating
#'   equation.
#' @returns A named list, empty when nothing forwarded describes the model.
#' @noRd
spec_ee_args <- function(ee_args) {
  ee_args[!names(ee_args) %in% per_observation_ee_args()]
}

#' Record what a formula fit was specified as
#'
#' Builds the list stored in a fitted estimator's `model_spec` property. Called
#' from `prepare_formula_psi()`, the one place that has seen the formula, the
#' model frame, the design, and the estimating equation together.
#'
#' The three fields that rebuild a design are `terms`, `xlevels`, and
#' `contrasts`, the same triple [stats::lm()] keeps and for the same reason: a
#' design matrix is not reproducible from data alone. A factor whose observed
#' levels are a subset of the fitted ones yields fewer columns, a non-default
#' contrast is not recoverable from the columns it produced, and a formula's
#' transformations have to be applied again. Given the triple,
#' `model.matrix(terms, model.frame(terms, data, xlev = xlevels), contrasts.arg
#' = contrasts)` reproduces the fitted design exactly. The terms carry the
#' response index and the offset attribute, so they also say which variable was
#' the response and which was offset by.
#'
#' The terms are read off the model frame rather than off the formula, which is
#' what makes a data-dependent transformation reproducible. `model.frame()`
#' records the fitted call of a term such as `poly(x, 2)` or `scale(x)` in the
#' terms' `predvars` attribute, coefficients and all, and evaluates that stored
#' call on any later data. Terms taken from the formula carry no such attribute,
#' so the basis would be recomputed on whatever rows it was handed and the
#' rebuilt design would silently disagree with the fitted one.
#'
#' `n_coef` is the number of design columns, which is not always the number of
#' parameters. `ee_glm()` under the gamma distribution solves for
#' `c(beta, log_shape)`, so a linear predictor formed as `X %*% theta` is
#' non-conformable and the coefficients are the leading `n_coef` entries. It is
#' counted from the design rather than taken from `appended_param_name()`, which
#' answers a different question: that table names the parameter an equation
#' appends, and returns `NULL` both for an equation that appends nothing and for
#' `ee_gformula()`, which puts its extra parameter first rather than last.
#' Counting design columns is exact wherever the table is silent.
#'
#' @param mf The model frame.
#' @param X The design matrix built from it.
#' @param y The response, after `coerce_formula_response()`.
#' @param .ee The estimating-equation function.
#' @param ee_args The evaluated `...` arguments forwarded to it, including the
#'   offset lifted out of the formula.
#' @param response_levels The levels of a factor or character response, or
#'   `NULL`.
#' @returns A named list.
#' @noRd
formula_model_spec <- function(mf, X, y, .ee, ee_args, response_levels) {
  model_terms <- stats::terms(mf)
  # `[[` on an absent name errors, and `$` partially matches, so the offset is
  # looked up by an exact membership test.
  offset <- if ("offset" %in% names(ee_args)) ee_args[["offset"]] else NULL

  list(
    terms = model_terms,
    xlevels = stats::.getXlevels(model_terms, mf),
    contrasts = attr(X, "contrasts"),
    ee = .ee,
    ee_spec_args = spec_ee_args(ee_args),
    n_coef = ncol(X),
    X = X,
    y = y,
    offset = offset,
    response_levels = response_levels
  )
}

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
# has not seen. Only the first kind may be applied to new data, so
# `spec_ee_args()` keeps those in `ee_spec_args` and `spec_obs_args()` puts the
# rest in `ee_obs_args`.
#
# Both halves are recorded, in separate fields, because a quantity of the fitted
# sample is sometimes needed to predict on data the fit has not seen.
# `plogit_predict()` is the case that settled it: the time grid a pooled logistic
# fit is defined on comes from the observed times and event indicators of the
# fitted sample, so predicting new covariate patterns on that grid needs the
# fitted `event` alongside the new design. Keeping the two halves apart is what
# stops such a value from being applied to new data as though it described the
# model. `offset` also has a field of its own, because `predict()` asks a
# question about it that no other per-observation argument raises: an offset in
# the terms is evaluated on `newdata` while one supplied through `...` cannot be,
# and that comparison is made before anything looks at the rest.
#
# Both halves hold whatever the fit happened to forward, which for most
# equations is a subset of what they take, so a reader of either has to treat
# absence as the ordinary case rather than as a fault. Every such read is
# written `args[["name"]]`, which is `NULL` for a name the fit did not forward,
# rather than `args$name`, which would partially match one recorded argument to
# another's name.
#
# Reading by name is also what bounds the record: the two halves partition the
# arguments that carry a name, and an argument forwarded without one is recorded
# in neither. `named_ee_args()` drops it before the split, for the reason given
# there.

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
#' An argument not on the list is treated as model specification, so this list is
#' the whole of the classification and neither direction of a mistake in it
#' announces itself at runtime. `predict()` never calls the estimating equation
#' again, and it reads either half by explicit name, so an argument wrongly kept
#' here is simply never asked for, and one wrongly dropped reads back as `NULL`,
#' which every reader takes for an argument the fit did not forward: a default
#' link, no appended parameter, or a report about a design that does not account
#' for the fit's parameters. Nothing says the specification was misfiled.
#'
#' That is why the list is pinned from both ends rather than only documented.
#' tests/testthat/test-model-spec-registry.R holds a second table classifying
#' every argument of every exported `ee_*()` function, derived from their
#' signatures, and checks it against this list and against the two functions that
#' read it. A new argument fails that file until it is classified, and a
#' classification this list disagrees with fails it too.
#'
#' `X` is the one observation-indexed design argument left off, where `W`, `V`,
#' `X1`, and `X0` are all on. The formula interface fills `X` by name, so an `X`
#' forwarded through `...` makes the call match that argument twice and the fit
#' fails before anything is classified. Nothing this list could say about it is
#' reachable, which the same test file pins.
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

#' The forwarded arguments that carry a name
#'
#' What both halves of the split are taken from, so that between them they
#' partition the named arguments and record nothing else.
#'
#' An argument forwarded without a name still reaches the estimating equation,
#' where `do.call()` matches it by position, but it cannot be recorded: both
#' halves are read by explicit name, so a nameless entry has no key any later
#' reader could ask for and keeping one would record a specification nothing can
#' reach.
#'
#' A name goes missing in three shapes, and this is where all of them are
#' answered the same way. `rlang::enquos()` gives a dot passed by position an
#' empty name rather than leaving it out, an argument list that named nothing has
#' no names vector at all, and a names vector assembled by hand can carry `NA`
#' where a name should be. `formula_ee_dots_names()` reads the first two as no
#' name supplied and so does this; the third is read the same way because a name
#' that is not there and a name that is missing are equally unusable as a key.
#' `nzchar()` alone would keep it, since `nzchar(NA)` is `TRUE`, and the halves
#' below split on `%in%`, which is `FALSE` for `NA` on both sides, so such an
#' entry would land in the model-specification half by falling through rather
#' than by belonging there. Leaving any of the three to the subscript would be an
#' accident rather than a rule: `!NULL %in% x` is `logical(0)`, and a list
#' subscripted by `logical(0)` is empty whatever it held.
#'
#' @param ee_args The evaluated `...` arguments forwarded to the estimating
#'   equation.
#' @returns The elements of `ee_args` that carry a name, in order.
#' @noRd
named_ee_args <- function(ee_args) {
  arg_names <- names(ee_args)
  if (is.null(arg_names)) {
    return(ee_args[0])
  }
  ee_args[!is.na(arg_names) & nzchar(arg_names)]
}

#' The model-specification subset of the arguments forwarded to an equation
#'
#' Drops the per-observation arguments named by `per_observation_ee_args()` and
#' keeps the rest of the named ones, which describe the model rather than the
#' sample it was fitted to.
#'
#' @param ee_args The evaluated `...` arguments forwarded to the estimating
#'   equation.
#' @returns A named list, empty when nothing forwarded describes the model.
#' @noRd
spec_ee_args <- function(ee_args) {
  named <- named_ee_args(ee_args)
  named[!names(named) %in% per_observation_ee_args()]
}

#' The per-observation subset of the arguments forwarded to an equation
#'
#' The complement of `spec_ee_args()` within the named arguments: those named by
#' `per_observation_ee_args()`, which carry one value per row of the data the
#' fit was made on.
#'
#' @param ee_args The evaluated `...` arguments forwarded to the estimating
#'   equation.
#' @returns A named list, empty when nothing forwarded varies by observation.
#' @noRd
spec_obs_args <- function(ee_args) {
  named <- named_ee_args(ee_args)
  named[names(named) %in% per_observation_ee_args()]
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
#' The model frame is recorded whole, beside the design built from it, because
#' the two answer different questions. The design is the numeric matrix the
#' estimating equation was handed, with factors already coded and
#' transformations already applied; the model frame is the variables the formula
#' named, which is what [stats::model.frame()] promises and what a row-wise
#' report of a fit is written against. It costs nothing to keep: the estimating
#' function closure is built in this frame and its environment holds the model
#' frame already, so the field is a second reference to an object the fit
#' carries either way rather than a second copy of it.
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

  list(
    terms = model_terms,
    xlevels = stats::.getXlevels(model_terms, mf),
    contrasts = attr(X, "contrasts"),
    ee = .ee,
    ee_spec_args = spec_ee_args(ee_args),
    ee_obs_args = spec_obs_args(ee_args),
    n_coef = ncol(X),
    model_frame = mf,
    X = X,
    y = y,
    offset = ee_args[["offset"]],
    response_levels = response_levels
  )
}

#' Rebuild the design matrix of data from the recorded specification
#'
#' The one place a design is built from the recorded triple, shared by
#' [`predict()`][deli-predict] and `model.matrix()` so that the design a
#' prediction is formed on and the design a fit reports are the same matrix.
#' Three choices make it worth writing once.
#'
#' The response is deleted from the terms, so covariate values that carry none
#' are enough. A design matrix names no response, and neither does the `newdata`
#' of a prediction have to.
#'
#' The recorded factor levels and contrasts are passed rather than the session's.
#' A factor observed at a subset of its fitted levels would otherwise yield fewer
#' columns, and a coding the fit recorded would be replaced by whatever
#' `getOption("contrasts")` says when the call is made, answering for a model that
#' was never fitted and giving no sign of it.
#'
#' Rows with missing values pass through under `na.pass` rather than being
#' dropped, so the design lines up with the rows of `data`.
#'
#' The model frame is returned beside the matrix because `predict()` reads the
#' offset off it, and building it twice would evaluate every term twice.
#'
#' @param spec The fit's `model_spec`.
#' @param data A data frame of covariate values.
#' @returns A list with `model_frame` and `X`.
#' @noRd
rebuild_design <- function(spec, data) {
  predictors <- stats::delete.response(spec$terms)
  mf <- stats::model.frame(
    predictors,
    data,
    na.action = stats::na.pass,
    xlev = spec$xlevels
  )

  list(
    model_frame = mf,
    X = stats::model.matrix(predictors, mf, contrasts.arg = spec$contrasts)
  )
}

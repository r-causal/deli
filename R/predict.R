# ---- Which parameters predict() may treat as coefficients --------------------
# `predict()` forms a linear predictor as `X %*% beta`, so it has to know which
# entries of `theta` are the coefficients on the recorded design and where they
# sit. `model_spec$n_coef` is the number of design columns, and on its own it
# answers neither question.
#
# Two built-in equations show why. Under one plan `ee_gformula()` orders its
# parameters `causal_mean, X_1, ..., X_k`, one more than the design has columns
# and so the same count as `ee_tobit()`, whose extra parameter trails instead.
# Slicing to `n_coef` would put the causal mean where the intercept belongs and
# shift every other coefficient by one, and nothing about the result would look
# wrong.
# `ee_additive_regression()` expands its design internally, so a fit with two
# design columns and one spline specification solves for four coefficients on
# the expanded basis. There `n_coef` is not a coefficient count at all, and
# slicing to it drops real coefficients rather than misplacing them.
#
# `appended_param_name()` cannot separate those cases, because it returns `NULL`
# both for an equation that appends nothing and for an equation it does not
# list. A `NULL` from it means "appends nothing" only for an equation already
# known to be in its table.
#
# `predict_link_name()` below is the table that makes `n_coef` sufficient. Every
# equation it lists puts one coefficient per design column at the front of
# `theta` and appends at most the single parameter `appended_param_name()`
# names, so the coefficients are `theta[seq_len(n_coef)]` and
# `predict_coef_index()` confirms the parameter count agrees before slicing. An
# equation the table does not list gets an error naming the reason, never a
# slice. That direction is deliberate: refusing to predict is visible, and a
# wrong slice is not.
#
# The table carries the second thing `predict()` needs as well, the inverse link
# that takes the linear predictor to the mean of the response. An equation whose
# linear predictor is not a link-scale mean is left out for that reason alone.
# `ee_aft()` puts its linear predictor on the log-time scale and `ee_tobit()`
# puts it on a latent scale, so `type = "response"` would not be the conditional
# mean of the response for either. An `ee_aft()` fit is still predictable, on
# the other surface this method carries: `times` asks for a survival measure at
# a set of times, which is the question that model answers. See
# R/predict-survival.R for the two surfaces and for why no fit is on both.
#
# ---- Specification arguments, not a replay -----------------------------------
# `predict()` reads the model from `model_spec$ee_spec_args` and never calls the
# estimating equation again. That split is what makes `newdata` possible:
# `model`, `distribution`, and `link` describe the model and apply to any data,
# while `weights`, `offset`, and `event` are one value per fitted observation
# and have no counterpart in data the model has not seen. An offset supplied
# through `...` therefore cannot travel onto `newdata`, and `predict()` aborts
# rather than predicting at an offset of zero. An offset written into the
# formula with `offset()` is part of the terms, so it is evaluated on `newdata`
# like any other variable.
#
# ---- Two argument names from base R ------------------------------------------
# `se.fit` and `level` are spelled as `stats::predict.lm()` spells them, rather
# than in this package's own vocabulary of snake_case names and `alpha`. This is
# a method on a base generic, and the shapes it returns are the shapes
# `predict.lm()` returns, so a reader who knows one knows the other. `confint()`
# takes `level` here for the same reason.

#' Predictions from a fitted deli estimator
#'
#' A method for [stats::predict()] that predicts from a fit made through the
#' formula interface, with sandwich standard errors and Wald confidence
#' intervals. It forms the linear predictor of a regression fit on either the
#' link or the response scale, and evaluates a survival measure at a set of
#' `times` for a fit of [ee_aft()] or [ee_plogit()].
#'
#' @details
#' # Two surfaces
#'
#' `times` chooses what is predicted. Without it, `predict()` returns the linear
#' predictor on the scale `type` names. With it, `predict()` returns the survival
#' measure `measure` names at each of the times, one prediction per row of the
#' design at each time.
#'
#' The two cover disjoint sets of fits. An equation reaches the first surface by
#' having a linear predictor that is a conditional mean of the response, and the
#' two equations on the second surface are there because theirs is not:
#' [ee_aft()] puts its linear predictor on the log-time scale, and [ee_plogit()]
#' has one linear predictor per person and time interval rather than one per
#' person. Supplying `type` beside `times`, or `measure`, `deriv_method`, or `dx`
#' without `times`, is an error rather than an argument silently ignored, since
#' no fit takes both sets.
#'
#' # The linear predictor
#'
#' The design for `newdata` is rebuilt through the terms, factor levels, and
#' contrasts the fit recorded, so a factor whose new values cover only some of
#' the fitted levels still produces the fitted set of columns, and a
#' data-dependent term such as `poly(x, 2)` or `scale(x)` is evaluated with the
#' coefficients it was fitted with rather than refitted to `newdata`. A factor
#' level the fit never saw, or a predictor `newdata` does not carry, is an
#' error.
#'
#' The standard error is the delta-method standard error of the linear
#' predictor, \eqn{\sqrt{\mathrm{diag}(X \hat{V} X^{T})}}, formed from the
#' coefficient block of the sandwich variance. On the response scale it is
#' scaled by the derivative of the inverse link, which is exact because the
#' inverse link is applied elementwise. The intervals are Wald intervals on the
#' scale asked for, so they are symmetric about `fit` on that scale rather than
#' transformed from the link scale.
#'
#' `predict()` needs to know which parameters are coefficients on the design and
#' what takes the linear predictor to the mean of the response. It supports
#' [ee_regression()], [ee_glm()], [ee_robust_regression()],
#' [ee_beta_regression()], and the five penalized regressions
#' [ee_bridge_regression()], [ee_ridge_regression()], [ee_lasso_regression()],
#' [ee_dlasso_regression()], and [ee_elasticnet_regression()], whose parameters
#' are one coefficient per design column followed by at most one parameter of
#' the outcome distribution. Any other estimating equation, and any fit built
#' from a `stacked_equations` function, is an error naming the reason;
#' [regression_predictions()] takes a design, estimates, and a covariance matrix
#' directly and can be used wherever this method declines.
#'
#' # Survival measures at a set of times
#'
#' `times` is supported for a fit of [ee_aft()] or of [ee_plogit()], and
#' `measure` names any of the measures [convert_survival_measures()] defines. The
#' point estimates are those of [aft_predictions_individual()] and
#' [plogit_predict()] for the same fit.
#'
#' [ee_survival_model()] has no surface here because the formula interface cannot
#' drive it: it takes no design matrix, and the interface always passes the one
#' it built. Predict from such a fit with [survival_predictions()].
#'
#' Each covariate pattern has its own variance, so each row gets its own
#' interval. The variance of the measure for a row at a time is the delta-method
#' variance \eqn{G \hat{V} G^{T}} of that one prediction, where \eqn{G} is its
#' row of the Jacobian of the whole grid with respect to the parameters, built by
#' `deriv_method`. A pooled logistic fit is predicted through
#' [plogit_predict()], whose matrix products cannot carry tangents, so
#' `deriv_method = "exact"` is available for an AFT fit alone. Asking for it on a
#' pooled logistic fit is an error whether or not a standard error is wanted,
#' since the Jacobian that could not be built is only built when one is.
#'
#' @param object A fitted `MEstimator` or `GMMEstimator` object made with the
#'   formula interface (after calling [estimate()]).
#' @param newdata A data frame of covariate values to predict at, or `NULL`
#'   (default) to predict the rows the fit was made on. An offset written into
#'   the formula with `offset()` is evaluated on `newdata`; an offset supplied
#'   through `...` at fit time is one value per fitted observation, so
#'   predicting on `newdata` is an error rather than a prediction at an offset
#'   of zero.
#' @param type Character string. `"link"` (default) returns the linear
#'   predictor; `"response"` applies the inverse link, giving the conditional
#'   mean of the response. Cannot be supplied beside `times`.
#' @param se.fit Logical. Return standard errors beside the predictions?
#'   Default `FALSE`.
#' @param interval Character string. `"none"` (default) or `"confidence"` for
#'   Wald intervals at `level`.
#' @param level The confidence level for `interval`. Default `0.95`. The
#'   critical value comes from the standard normal distribution, or from the
#'   t-distribution with \eqn{n - p} degrees of freedom when a
#'   `finite_correction` is set on the fit, matching
#'   [`confint()`][deli-generics].
#' @param times Numeric vector of times to predict a survival measure at, or
#'   `NULL` (default) to predict the linear predictor instead. Supported for a
#'   fit of [ee_aft()] or [ee_plogit()].
#' @param measure Character string naming the survival measure, one of
#'   `"survival"` (default), `"risk"`, `"cumulative_hazard"`, `"hazard"`, or
#'   `"density"`; see [convert_survival_measures()]. Only meaningful beside
#'   `times`.
#' @param deriv_method Character string for the derivative method used to build
#'   the Jacobian of the survival measure. One of `"capprox"` (central
#'   difference, the default), `"fapprox"`, `"bapprox"`, or `"exact"`
#'   (forward-mode automatic differentiation, available for an AFT fit and an
#'   error for a pooled logistic one). Only meaningful beside `times`.
#' @param dx Numeric step size for the finite-difference methods, ignored when
#'   `deriv_method = "exact"`. Default `1e-9`. Only meaningful beside `times`.
#'   Must be a single positive finite number, which is checked wherever it is
#'   supplied beside `times`, including a prediction asking for no standard
#'   error, which takes no step at all.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns
#' Without `times`, a named numeric vector of predictions, one per row of
#' `newdata` or of the fitted design, and with `interval = "confidence"` a matrix
#' with columns `"fit"`, `"lwr"`, and `"upr"`.
#'
#' With `times`, a data frame with one row per row of the design and time, every
#' time for the first row before any time for the second, and columns `.row`, the
#' row label of the design, `time`, and `fit`. With `interval = "confidence"` it
#' also has `"lwr"` and `"upr"`.
#'
#' With `se.fit = TRUE`, either shape becomes a list whose `fit` element is
#' whichever of the two the other arguments call for and whose `se.fit` element
#' is a numeric vector of standard errors.
#'
#' @seealso [deli-augment], which returns the linear-predictor predictions as
#'   columns beside the data; [regression_predictions()], which computes the same
#'   quantities from a design matrix, a vector of estimates, and a covariance
#'   matrix, for fits this method does not cover; and
#'   [aft_predictions_individual()], [aft_predictions_function()],
#'   [plogit_predict()], and [survival_predictions()], which compute the survival
#'   measures from estimates directly.
#'
#' @examples
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' head(predict(fit))
#'
#' # New covariate patterns, with sandwich standard errors.
#' at <- data.frame(wt = c(2, 3, 4), hp = 110)
#'
#' predict(fit, newdata = at, se.fit = TRUE)
#'
#' predict(fit, newdata = at, interval = "confidence")
#'
#' # A logistic fit predicts log-odds by default and probabilities on the
#' # response scale.
#' vs_fit <- m_estimate(vs ~ wt, data = mtcars, .ee = ee_regression,
#'                      model = "logistic")
#'
#' predict(vs_fit, newdata = data.frame(wt = c(2, 3)), type = "response")
#'
#' # A pooled logistic fit has no linear predictor to put on a scale, and
#' # predicts a survival measure at a set of times instead. Every row of the
#' # design gets its own interval.
#' bladder <- collett_bladder
#' bladder$novel <- bladder$treat - 1
#' k <- length(unique(bladder$time[bladder$delta == 1]))
#'
#' plogit_fit <- m_estimate(
#'   time ~ novel + init + size - 1, data = bladder, .ee = ee_plogit,
#'   event = delta, init = c(rep(0, 3), -4, rep(0, k - 1))
#' )
#'
#' head(predict(plogit_fit, times = c(12, 24), interval = "confidence"))
#'
#' @name deli-predict
#' @importFrom stats predict
NULL

# The S3 class of an S7 object is package-qualified, so the method below is
# bound to a backticked name and registered through a NAMESPACE S3method
# directive; see the comment above the accessors in R/generics-s3.R.

# ---- predict -----------------------------------------------------------------

#' @rdname deli-predict
#' @export
`predict.deli::deli_estimator` <- function(
  object,
  newdata = NULL,
  type = c("link", "response"),
  se.fit = FALSE,
  interval = c("none", "confidence"),
  level = 0.95,
  times = NULL,
  measure = "survival",
  deriv_method = "capprox",
  dx = 1e-9,
  ...
) {
  # Dispatch leaves the generic's frame directly beneath this method, so
  # `sys.call(-1)` is the `predict()` call the caller wrote and the error names
  # it. The default would name this method's own frame, spelled
  # `predict.deli::deli_estimator()`, and `rlang::caller_env()` would name
  # whichever function called `predict()`, neither of which is where the
  # offending argument was written. Nothing beneath this frame gives `NULL`,
  # which reports the message with no call at all.
  rlang::check_dots_empty(call = sys.call(-1))
  # A prediction is formed from the estimates, so that is what is required here.
  # The variance is required where a standard error is actually asked for, below
  # and in `predict_survival()`, so a fit whose variance could not be built still
  # predicts and still says why it cannot put an interval around it.
  check_has_estimates(object)
  # Asked before `match.arg()`, which fills `type` in and so erases the
  # difference between a caller who wrote it and one who did not.
  check_predict_surface(
    times,
    type_given = !missing(type),
    survival_given = c(
      measure = !missing(measure),
      deriv_method = !missing(deriv_method),
      dx = !missing(dx)
    )
  )
  interval <- match.arg(interval)
  check_level(level)

  if (!is.null(times)) {
    check_predict_times(times)
    check_predict_measure(measure)
    # Validated here rather than where they are used, since they are only used
    # when a standard error is asked for and a misspelled method or an
    # impossible step is worth reporting either way.
    deriv_method <- check_deriv_method(deriv_method)
    check_dx(dx)
    return(predict_survival(
      object,
      newdata = newdata,
      times = as.numeric(times),
      measure = measure,
      se.fit = se.fit,
      interval = interval,
      level = level,
      deriv_method = deriv_method,
      dx = dx
    ))
  }

  type <- match.arg(type)
  link <- resolve_predict_link(object, "predict")
  spec <- object@model_spec

  index <- predict_coef_index(spec, object@n_params)
  beta <- object@theta[index]

  design <- predict_design(spec, newdata)
  X <- design$X

  eta <- as.numeric(X %*% beta)
  if (!is.null(design$offset)) {
    eta <- eta + as.numeric(design$offset)
  }

  if (type == "response") {
    inverse <- inverse_link(eta, link)
    fit <- inverse$mu
  } else {
    inverse <- NULL
    fit <- eta
  }
  names(fit) <- rownames(X)

  # Everything below reads the variance, and nothing above it does, so a fit
  # that has none answers this far and no further.
  if (!se.fit && interval == "none") {
    return(fit)
  }
  check_estimated(object)

  covariance <- object@variance[index, index, drop = FALSE]
  # The variance of the linear predictor is the diagonal of X V X'. It is
  # written the way regression_predictions() writes it, rather than in any other
  # algebraically equal form, so that the two agree to the last bit.
  eta_variance <- as.numeric(rowSums((X %*% covariance) * X))
  fit_variance <- if (is.null(inverse)) {
    eta_variance
  } else {
    # The inverse link is applied one observation at a time, so its Jacobian is
    # diagonal and the delta-method variance of the mean is the squared
    # derivative times the variance of the linear predictor, with no
    # approximation beyond the delta method itself.
    inverse$dmu^2 * eta_variance
  }

  se <- sqrt(fit_variance)
  names(se) <- rownames(X)

  if (interval == "confidence") {
    cv <- critical_value(object, 1 - level)
    fit <- cbind(fit = fit, lwr = fit - cv * se, upr = fit + cv * se)
  }

  if (!se.fit) {
    return(fit)
  }
  list(fit = fit, se.fit = se)
}

# ---- The supported equations -------------------------------------------------

#' The inverse link of an estimating equation predict() can predict from
#'
#' Names the link that takes an equation's linear predictor to the mean of its
#' response, for the equations whose parameters are one coefficient per design
#' column followed by at most the parameter `appended_param_name()` names.
#' `NULL` for every other equation, including every equation not listed, which
#' is what `predict()` refuses on.
#'
#' Membership carries the layout claim as well as the link, which is why the two
#' answers are not separated. An equation reaches this table only when both hold
#' of it, so a caller that has a link from here may slice the coefficients off
#' the front of `theta`.
#'
#' The three models [ee_regression()] takes are the three links they name:
#' `"linear"` is the identity link, `"logistic"` the logit, and `"poisson"` the
#' log. Returning a link for those rather than the model itself lets every scale
#' run through `inverse_link()`, which supplies the derivative of the inverse
#' link as well as its value, where `model_transform()` supplies only the value.
#' `test-predict.R` pins the two against each other so the mapping cannot drift.
#'
#' @param .ee The estimating-equation function the fit was made with.
#' @param ee_args The model-specification arguments recorded for the fit.
#' @returns A string naming a link `inverse_link()` accepts, or `NULL`.
#' @noRd
predict_link_name <- function(.ee, ee_args) {
  if (!is.function(.ee)) {
    return(NULL)
  }
  if (identical(.ee, ee_glm)) {
    return(link_arg_name(ee_args[["link"]]))
  }
  if (identical(.ee, ee_beta_regression)) {
    # The mean model is logit by construction; the equation takes no link.
    return("logit")
  }
  model_equations <- list(
    ee_regression,
    ee_robust_regression,
    ee_ridge_regression,
    ee_bridge_regression,
    ee_lasso_regression,
    ee_dlasso_regression,
    ee_elasticnet_regression
  )
  if (any(vapply(model_equations, identical, logical(1), y = .ee))) {
    return(model_link_name(ee_args[["model"]]))
  }
  NULL
}

#' The link named by a `model` argument
#'
#' @param model The `model` argument recorded for the fit.
#' @returns A string, or `NULL` when the argument names no supported model.
#' @noRd
model_link_name <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    return(NULL)
  }
  switch(
    tolower(model),
    linear = "identity",
    logistic = "logit",
    poisson = "log"
  )
}

#' The link named by a `link` argument
#'
#' Passed on unchanged for `inverse_link()` to accept or reject, since the set
#' of links it takes is its own to state.
#'
#' @param link The `link` argument recorded for the fit.
#' @returns A string, or `NULL` when the argument is not one.
#' @noRd
link_arg_name <- function(link) {
  if (!is.character(link) || length(link) != 1L || is.na(link)) {
    return(NULL)
  }
  link
}

#' Where the design coefficients sit in a fit's parameter vector
#'
#' Returns the positions of the coefficients on the recorded design, having
#' first confirmed that the design and the parameter the estimating equation
#' appends account for every parameter the fit solved for. They will not when a
#' caller supplied an `init` longer or shorter than the equation's layout, so
#' the check is what stops a design of the wrong width from being multiplied
#' into a slice of the wrong parameters.
#'
#' Only equations `appended_param_name()` lists reach here, whether through
#' `predict_link_name()` or through the AFT branch of `predict_survival()`, so a
#' `NULL` from it means the equation appends nothing rather than that the
#' equation is unknown to it.
#'
#' @param spec The fit's `model_spec`.
#' @param n_params The number of parameters the fit solved for.
#' @returns An integer vector of positions in `theta`.
#' @noRd
predict_coef_index <- function(spec, n_params) {
  n_coef <- spec$n_coef
  appended <- appended_param_name(spec$ee, spec$ee_spec_args)
  expected <- if (is.null(appended)) n_coef else n_coef + 1L

  if (!identical(as.integer(n_params), as.integer(expected))) {
    cli::cli_abort(
      c(
        "The recorded design does not account for every parameter of this fit.",
        "i" = "The model matrix has {n_coef} column{?s} and the fit solved for
               {n_params} parameter{?s}, where {expected} {?was/were} expected.",
        "i" = "{.fn predict} cannot tell which parameters are the coefficients,
               so it makes no prediction rather than a wrong one."
      ),
      call = NULL
    )
  }
  seq_len(n_coef)
}

# ---- The design to predict on ------------------------------------------------

#' The design and offset to predict at
#'
#' Returns the fitted design when `newdata` is absent, and otherwise the design
#' `rebuild_design()` builds from the terms, factor levels, and contrasts the fit
#' recorded, the same triple [stats::lm()] keeps and the same way
#' [stats::predict.lm()] uses it. The offset is the question this adds to that
#' one, and it is why the model frame comes back beside the matrix.
#'
#' @param spec The fit's `model_spec`.
#' @param newdata A data frame, or `NULL`.
#' @returns A list with `X` and `offset`.
#' @noRd
predict_design <- function(spec, newdata) {
  if (is.null(newdata)) {
    return(list(X = spec$X, offset = spec$offset))
  }

  design <- rebuild_design(spec, newdata)

  # An offset in the terms is evaluated on newdata like any other variable. One
  # that is not in the terms reached the fit through `...`, and is one value per
  # fitted observation rather than a rule newdata can be put through.
  offset <- stats::model.offset(design$model_frame)
  if (is.null(offset) && !is.null(spec$offset)) {
    abort_predict_dots_offset()
  }

  list(X = design$X, offset = offset)
}

# ---- Validation and conditions -----------------------------------------------

#' The link a fit can be predicted through, or an error naming why it cannot
#'
#' Answers the two questions that decide whether a fit can be predicted from at
#' all: whether it records a model specification, and whether its estimating
#' equation is one whose linear predictor `predict_link_name()` can name a link
#' for. Both are asked here, in one place, so that [`predict()`][deli-predict],
#' `fitted()`, `residuals()`, and `augment()` refuse the same fits on the same
#' terms.
#'
#' `fn` is what makes it worth sharing. Every one of those functions reaches the
#' supported-equation table through this helper, and a user who called
#' `augment()` should be told that `augment()` declined rather than be handed a
#' message about a function they did not call.
#'
#' @param object A fitted estimator.
#' @param fn The name of the function the user called.
#' @returns A string naming a link `inverse_link()` accepts.
#' @noRd
resolve_predict_link <- function(object, fn = "predict") {
  spec <- object@model_spec
  if (is.null(spec)) {
    abort_predict_no_model_spec(fn)
  }
  link <- predict_link_name(spec$ee, spec$ee_spec_args)
  if (is.null(link)) {
    abort_predict_unsupported_ee(spec$ee, fn)
  }
  link
}

#' Validate a confidence level
#'
#' @param level The level supplied.
#' @param arg The name the caller gave it. `augment()` spells it `conf.level`,
#'   and the message names the argument the user wrote rather than the one this
#'   package passes it on as.
#' @returns Invisible `NULL`. Raises an error if `level` is not a single number
#'   strictly between 0 and 1.
#' @noRd
check_level <- function(level, arg = "level") {
  if (
    !is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1
  ) {
    cli::cli_abort(
      "{.arg {arg}} must be a single number between 0 and 1 (exclusive).",
      call = NULL
    )
  }
  invisible(NULL)
}

#' @noRd
abort_predict_no_model_spec <- function(fn = "predict") {
  cli::cli_abort(
    c(
      "{.fn {fn}} needs a fit made through the formula interface.",
      "i" = "This fit was built from a {.arg stacked_equations} function, which
             records no formula, design, or link to predict from.",
      "i" = "Compute predictions from the estimates directly with
             {.fn regression_predictions}."
    ),
    call = NULL
  )
}

#' @noRd
abort_predict_unsupported_ee <- function(.ee, fn = "predict") {
  if (identical(.ee, ee_gformula)) {
    cli::cli_abort(
      c(
        "{.fn {fn}} does not support a fit of {.fn ee_gformula}.",
        "i" = "Its parameters are the causal mean or means it estimates
               followed by the outcome-model coefficients, so the leading
               parameters are not coefficients on the design.",
        "i" = "Predict from the outcome model with
               {.fn regression_predictions}, using the trailing coefficients
               and the matching block of {.fn vcov}."
      ),
      call = NULL
    )
  }
  if (identical(.ee, ee_aft)) {
    cli::cli_abort(
      c(
        "{.fn {fn}} does not form a linear predictor for a fit of
         {.fn ee_aft}.",
        "i" = "An accelerated failure time model puts its linear predictor on
               the log-time scale rather than on a link scale, so
               {.code type = \"response\"} would not be the conditional mean of
               the response.",
        "i" = "Predict a survival measure at a set of times instead, with
               {.code predict(object, times = , measure = )}."
      ),
      call = NULL
    )
  }
  if (identical(.ee, ee_plogit)) {
    cli::cli_abort(
      c(
        "{.fn {fn}} does not form a linear predictor for a fit of
         {.fn ee_plogit}.",
        "i" = "A pooled logistic model has one linear predictor per person and
               time interval rather than one per person, so there is no
               prediction to report against the rows of the data.",
        "i" = "Predict a survival measure at a set of times instead, with
               {.code predict(object, times = , measure = )}."
      ),
      call = NULL
    )
  }
  if (identical(.ee, ee_additive_regression)) {
    cli::cli_abort(
      c(
        "{.fn {fn}} does not support a fit of {.fn ee_additive_regression}.",
        "i" = "Its coefficients sit on the spline basis
               {.fn additive_design_matrix} expands the design into, not on the
               columns of the model matrix the formula built.",
        "i" = "Expand the new data with {.fn additive_design_matrix} and predict
               from the result with {.fn regression_predictions}."
      ),
      call = NULL
    )
  }
  cli::cli_abort(
    c(
      "{.fn {fn}} does not support this fit's estimating equation.",
      "i" = "It supports {.fn ee_regression}, {.fn ee_glm},
             {.fn ee_robust_regression}, {.fn ee_beta_regression}, and the five
             penalized regressions {.fn ee_bridge_regression},
             {.fn ee_ridge_regression}, {.fn ee_lasso_regression},
             {.fn ee_dlasso_regression}, and {.fn ee_elasticnet_regression},
             whose parameters are one coefficient per design column followed by
             at most one parameter of the outcome distribution.",
      "i" = "For any other equation, build the design yourself and predict with
             {.fn regression_predictions}."
    ),
    call = NULL
  )
}

#' @noRd
abort_predict_dots_offset <- function() {
  cli::cli_abort(
    c(
      "The offset this fit was made with cannot be evaluated on
       {.arg newdata}.",
      "i" = "It was supplied through {.arg ...}, which is one value per fitted
             observation and has no counterpart in data the fit has not seen.",
      "i" = "Write the offset into the formula with {.fn offset} so that it is
             evaluated on {.arg newdata} like any other variable."
    ),
    call = NULL
  )
}

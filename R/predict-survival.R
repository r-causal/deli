# ---- Two surfaces on one generic ---------------------------------------------
# `predict()` answers two different questions for two disjoint sets of fits, and
# `times` is the switch between them.
#
# Without `times` it forms the linear predictor `X %*% beta` and puts it on the
# link or the response scale, which is what `type` names. That surface needs an
# equation whose linear predictor is a link-scale conditional mean, and
# `predict_link_name()` in R/predict.R is the table of those.
#
# With `times` it evaluates a survival measure at each of them, which is what
# `measure` names. That surface needs an equation that defines a survival
# function, and `predict_survival_kind()` below is the table of those.
#
# No equation is in both tables, and that is not an accident of which ones have
# been written down. An equation reaches the first table by having a linear
# predictor that is a conditional mean of the response, and the two equations in
# the second table are there precisely because theirs is not: `ee_aft()` puts its
# linear predictor on the log-time scale, and `ee_plogit()` has one value per
# person and time interval rather than one per person. So the surfaces partition
# the fits rather than overlapping on any of them, and a fit that accepts `times`
# is a fit that refuses `type`.
#
# That is why supplying both is an error rather than one of them winning. There
# is no fit for which both are meaningful, so a call carrying both was written
# under a mistaken belief about the fit, and the argument that would have been
# ignored is the one that says which belief. The check is on whether the caller
# wrote `type`, not on its value, because `type = "link"` is as much a statement
# about the scale as `type = "response"` is. `measure`, `deriv_method`, and `dx`
# are refused without `times` for the same reason in the other direction.
#
# ---- Which equations have a survival surface ---------------------------------
# `ee_aft()` and `ee_plogit()`. Two more are worth naming for the reason they are
# absent.
#
# `ee_survival_model()` cannot be fitted through the formula interface at all. It
# takes no design matrix, and `prepare_formula_psi()` always passes the design it
# built as `X`, so driving it from a formula is an unused-argument error before
# anything is estimated. There is therefore no fit of it to route, and
# `survival_predictions()` is how its predictions are made.
#
# `ee_mlogit()` is absent for the plainer reason that a multinomial model defines
# no survival function. It is reachable from the formula interface, which its
# indicator-matrix response might suggest it is not: `coerce_formula_response()`
# rejects only a factor or character response with other than two levels, and a
# matrix is neither, so `cbind(y1, y2, y3) ~ x` passes the matrix through into
# `y` and fits. Fits of it therefore exist, and `predict_survival_kind()` returns
# `NULL` for it rather than relying on there being none.
#
# ---- One transform, one Jacobian, one variance per row -----------------------
# `aft_predictions_function()` requires a single covariate pattern because each
# pattern has its own variance, and it forms the whole delta-method covariance
# across times for that one pattern and keeps only the diagonal. Doing that once
# per row would be n Jacobians for n patterns, and forming n covariance matrices
# to discard everything off their diagonals.
#
# The work is arranged the other way here. One transform returns the measure at
# every (row, time) pair, so one Jacobian covers the whole grid, and the variance
# of each pair is `rowSums((G %*% Sigma) * G)`, the diagonal of the delta-method
# covariance without the matrix it is the diagonal of. The Jacobian costs one
# evaluation of the transform per parameter, or two under a central difference,
# whatever n and k are, and the quadratic form costs one pass over the grid. The
# per-row variances are exactly the ones the per-row calls would produce: row
# `(i, j)` of the Jacobian is the derivative of that pair's measure and of
# nothing else, so no pattern's variance is borrowed by another.

#' The survival surface an estimating equation has, if it has one
#'
#' Names the family of survival predictions `predict()` can make from a fit of
#' an equation, and `NULL` for every equation with none, including every
#' equation not listed. Membership is by identity, so a user function that wraps
#' one of these equations is not routed: the wrapper may reorder or reinterpret
#' the parameters, and nothing about it says whether it did.
#'
#' @param .ee The estimating-equation function the fit was made with.
#' @returns `"aft"`, `"plogit"`, or `NULL`.
#' @noRd
predict_survival_kind <- function(.ee) {
  if (!is.function(.ee)) {
    return(NULL)
  }
  if (identical(.ee, ee_aft)) {
    return("aft")
  }
  if (identical(.ee, ee_plogit)) {
    return("plogit")
  }
  NULL
}

#' The survival surface a fit can be predicted on, or an error naming why not
#'
#' The counterpart of `resolve_predict_link()` for the `times` surface, asking
#' the same two questions: whether the fit records a model specification, and
#' whether its estimating equation is one this surface covers.
#'
#' @param object A fitted estimator.
#' @param fn The name of the function the user called.
#' @returns `"aft"` or `"plogit"`.
#' @noRd
resolve_predict_survival_kind <- function(object, fn = "predict") {
  spec <- object@model_spec
  if (is.null(spec)) {
    abort_predict_no_model_spec(fn)
  }
  kind <- predict_survival_kind(spec$ee)
  if (is.null(kind)) {
    abort_predict_no_survival_measure(fn)
  }
  kind
}

# ---- The survival surface ----------------------------------------------------

#' Predictions of a survival measure at a set of times
#'
#' The body of [`predict()`][deli-predict] when `times` is supplied. Everything
#' the two surfaces share is decided by the caller and passed in.
#'
#' @param object A fitted estimator.
#' @param newdata A data frame, or `NULL`.
#' @param times The times to predict at.
#' @param measure The survival measure to return.
#' @param se.fit Return standard errors beside the predictions?
#' @param interval `"none"` or `"confidence"`, already matched.
#' @param level The confidence level.
#' @param deriv_method The derivative method, already normalized.
#' @param dx The finite-difference step.
#' @returns A data frame, or a list of one and a vector of standard errors.
#' @noRd
predict_survival <- function(
  object,
  newdata,
  times,
  measure,
  se.fit,
  interval,
  level,
  deriv_method,
  dx
) {
  kind <- resolve_predict_survival_kind(object, "predict")
  spec <- object@model_spec
  if (kind == "plogit") {
    # Both refusals are made here, before the design is rebuilt and before a
    # standard error is asked for, so that neither depends on how the rest of
    # the call was written.
    #
    # An offset asked about now is the one the fit recorded, whichever way it
    # was supplied. Left to `predict_design()`, a fit whose offset came through
    # `...` would be told to write it into the formula instead, and a fit that
    # had done so would be refused here all the same.
    if (!is.null(spec$offset)) {
      abort_predict_plogit_offset()
    }
    # The Jacobian is only built when a standard error is wanted, so a call
    # without one would otherwise accept a derivative method that no call with
    # one can use.
    if (deriv_method == "exact") {
      abort_predict_plogit_exact()
    }
  }
  design <- predict_design(spec, newdata)

  built <- if (kind == "aft") {
    # The same layout check the linear-predictor surface makes: the coefficients
    # are the leading parameters and at most the one parameter of the error
    # distribution follows them, so a fit that solved for some other number is
    # refused rather than sliced. The scale is read off the last parameter, as
    # `aft_predictions_individual()` reads it.
    aft_prediction_transform(
      X = design$X,
      offset = design$offset,
      times = times,
      distribution = spec$ee_spec_args[["distribution"]],
      measure = measure,
      beta_index = predict_coef_index(spec, object@n_params),
      scale_index = object@n_params
    )
  } else {
    plogit_prediction_transform(
      spec = spec,
      X = design$X,
      times = times,
      measure = measure
    )
  }

  rows <- predict_row_labels(design$X)
  fit <- as.numeric(built$g(object@theta))[built$order]
  frame <- data.frame(
    .row = rep(rows, each = length(times)),
    time = rep(times, times = length(rows)),
    fit = fit,
    stringsAsFactors = FALSE
  )

  if (!se.fit && interval == "none") {
    return(frame)
  }

  # One call differentiates the transform and so evaluates it once or twice per
  # parameter; the scope collapses the repeats of any warning it raises, as
  # `delta_method()` does for the same reason.
  jacobian <- without_repeated_warnings(
    transform_jacobian(object@theta, built$g, deriv_method, dx)
  )
  # The diagonal of the delta-method covariance, written the way
  # `regression_predictions()` writes its quadratic form. pmax() at zero guards
  # the finite-difference paths, where step-size cancellation can drive a
  # near-zero variance slightly negative; it is a no-op under exact autodiff.
  variance <- rowSums((jacobian %*% object@variance) * jacobian)
  se <- sqrt(pmax(variance[built$order], 0))

  if (interval == "confidence") {
    cv <- critical_value(object, 1 - level)
    frame$lwr <- frame$fit - cv * se
    frame$upr <- frame$fit + cv * se
  }

  if (!se.fit) {
    return(frame)
  }
  list(fit = frame, se.fit = se)
}

#' The row labels the predictions are reported against
#'
#' @param X The design predicted on.
#' @returns A character vector, one entry per row.
#' @noRd
predict_row_labels <- function(X) {
  labels <- rownames(X)
  if (is.null(labels)) {
    return(as.character(seq_len(nrow(X))))
  }
  as.character(labels)
}

#' The order that takes a transform's output to the reported order
#'
#' The frame reports one row per observation and time, every time for the first
#' observation before any time for the second. A transform that produces the
#' grid a whole time at a time is in the other order, so its output is indexed
#' by this before being reported.
#'
#' @param n The number of observations.
#' @param k The number of times.
#' @returns An integer vector of length `n * k`.
#' @noRd
observation_major_order <- function(n, k) {
  as.integer(t(matrix(seq_len(n * k), nrow = n, ncol = k)))
}

#' The AFT measure across a design and a set of times, as one transform
#'
#' Returns the transform whose Jacobian gives the variance of every prediction,
#' together with the order that takes its output to the reported one.
#'
#' The linear predictor is accumulated one design column at a time rather than
#' through a matrix product, and that is what lets a single transform serve both
#' derivative methods. Under exact differentiation `th` is a tangent-carrying
#' pair vector, and reading a matrix product back out needs `as.numeric()`, which
#' a tangent-carrying value refuses; multiplying a column of the design by one
#' element of `th` is arithmetic, so the tangents survive it. The column is a
#' whole column, so the accumulation stays vectorized over observations and the
#' transform costs one pass per time rather than one per observation.
#'
#' @param X The design to predict on.
#' @param offset The offset for those rows, or `NULL`.
#' @param times The times to predict at.
#' @param distribution The error distribution recorded for the fit.
#' @param measure The measure to return.
#' @param beta_index The positions in `theta` of the coefficients on `X`, as
#'   `predict_coef_index()` returns them, so that the design is multiplied by
#'   exactly the parameters that check identified.
#' @param scale_index The position of the log inverse scale, which is the last
#'   parameter for every distribution but the exponential, where the scale is
#'   fixed at one and nothing is read from it.
#' @returns A list with `g` and `order`.
#' @noRd
aft_prediction_transform <- function(
  X,
  offset,
  times,
  distribution,
  measure,
  beta_index,
  scale_index
) {
  distribution <- check_predict_distribution(distribution)

  g <- function(th) {
    xbeta <- 0
    for (j in seq_along(beta_index)) {
      xbeta <- xbeta + X[, j] * th[beta_index[j]]
    }
    if (!is.null(offset)) {
      xbeta <- xbeta + as.numeric(offset)
    }
    sigma <- if (distribution == "exponential") 1 else exp(-th[scale_index])
    # c() rather than vapply(): under exact autodiff the per-time values are
    # tangent-carrying containers, not plain doubles, so no numeric template
    # can hold them.
    do.call(
      c,
      lapply(times, function(t) {
        aft_measure_at_time(t, xbeta, sigma, distribution, measure)
      })
    )
  }

  list(g = g, order = observation_major_order(nrow(X), length(times)))
}

#' The pooled logistic measure across a design and a set of times
#'
#' [plogit_predict()] already answers for the whole grid, so the transform is a
#' call to it. The times and event indicators are the fitted ones even when the
#' design is not: the time grid a pooled logistic fit is defined on is a property
#' of the sample it was fitted to, and a new covariate pattern is predicted on
#' that grid rather than on one of its own.
#'
#' Its return is times-by-observations, whose column-major flattening is already
#' one observation at a time, so nothing is reordered.
#'
#' There is no `offset` argument because [plogit_predict()] has none, and a
#' prediction made through it would silently be the prediction at an offset of
#' zero. `predict_survival()` refuses such a fit before reaching here.
#'
#' @param spec The fit's `model_spec`.
#' @param X The design to predict on.
#' @param times The times to predict at.
#' @param measure The measure to return.
#' @returns A list with `g` and `order`.
#' @noRd
plogit_prediction_transform <- function(spec, X, times, measure) {
  fitted_time <- spec$y
  fitted_event <- spec$ee_obs_args[["event"]]
  time_design <- spec$ee_spec_args[["S"]]
  unique_times <- spec$ee_spec_args[["unique_times"]]

  g <- function(th) {
    as.numeric(plogit_predict(
      th,
      time = fitted_time,
      event = fitted_event,
      X = X,
      S = time_design,
      times_to_predict = times,
      measure = measure,
      unique_times = unique_times
    ))
  }

  list(g = g, order = seq_len(nrow(X) * length(times)))
}

# ---- Validation and conditions -----------------------------------------------

#' Validate the times a survival measure is asked for at
#'
#' @param times The times supplied.
#' @returns Invisible `NULL`. Raises an error unless `times` is a non-empty
#'   vector of finite numbers.
#' @noRd
check_predict_times <- function(times) {
  if (!is.numeric(times) || length(times) == 0L || !all(is.finite(times))) {
    cli::cli_abort(
      c(
        "{.arg times} must be a numeric vector of finite times, not
         {.obj_type_friendly {times}}.",
        "i" = "It names the times a survival measure is evaluated at, one
               prediction per row of the design at each of them."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

#' Validate a survival measure
#'
#' Only the shape is checked here. Which measures exist is
#' [convert_survival_measures()]'s to state, and it states it once for every
#' function in the package that takes one.
#'
#' @param measure The measure supplied.
#' @returns Invisible `NULL`. Raises an error unless `measure` is a single
#'   string.
#' @noRd
check_predict_measure <- function(measure) {
  if (!is.character(measure) || length(measure) != 1L || is.na(measure)) {
    cli::cli_abort(
      "{.arg measure} must be a single string, not
       {.obj_type_friendly {measure}}.",
      call = NULL
    )
  }
  invisible(NULL)
}

#' Validate the distribution recorded for an AFT fit
#'
#' @param distribution The `distribution` argument recorded for the fit.
#' @returns The distribution, lowercased.
#' @noRd
check_predict_distribution <- function(distribution) {
  if (
    !is.character(distribution) ||
      length(distribution) != 1L ||
      is.na(distribution)
  ) {
    cli::cli_abort(
      c(
        "This fit records no error distribution to predict with.",
        "i" = "{.fn ee_aft} takes the distribution as a single string, and the
               survival measure at a time is defined by it."
      ),
      call = NULL
    )
  }
  tolower(distribution)
}

#' Validate that the arguments supplied name one surface
#'
#' `times` selects the surface, so `type` belongs to a call without it and
#' `measure`, `deriv_method`, and `dx` to a call with it. Whether each was
#' written is what is asked, not what it was set to: a default written out is
#' still a statement about the call.
#'
#' @param times The `times` argument.
#' @param type_given Was `type` supplied?
#' @param survival_given A named logical vector saying which of the arguments
#'   that belong to the `times` surface were supplied.
#' @returns Invisible `NULL`. Raises an error if the two surfaces are mixed.
#' @noRd
check_predict_surface <- function(times, type_given, survival_given) {
  if (is.null(times)) {
    supplied <- names(survival_given)[survival_given]
    if (length(supplied) > 0L) {
      cli::cli_abort(
        c(
          "{.arg {supplied}} {?is/are} only meaningful beside {.arg times}.",
          "i" = "{.arg times} asks for a survival measure at each of a set of
                 times, which is the only prediction {.arg {supplied}}
                 describe{?s/}.",
          "i" = "Without it {.fn predict} returns the linear predictor, on the
                 scale {.arg type} names."
        ),
        call = NULL
      )
    }
    return(invisible(NULL))
  }
  if (type_given) {
    cli::cli_abort(
      c(
        "{.arg type} and {.arg times} cannot both be supplied.",
        "i" = "{.arg type} names a scale of the linear predictor, and a
               survival measure at a time is not on one.",
        "i" = "Drop {.arg type} for a survival measure, or {.arg times} for the
               linear predictor. No fit takes both: an equation has a survival
               measure only when its linear predictor is not a conditional
               mean."
      ),
      call = NULL
    )
  }
  invisible(NULL)
}

#' @noRd
abort_predict_no_survival_measure <- function(fn = "predict") {
  cli::cli_abort(
    c(
      "{.fn {fn}} has no survival measure for this fit's estimating equation.",
      "i" = "{.arg times} and {.arg measure} are supported for {.fn ee_aft} and
             {.fn ee_plogit}, the two equations the formula interface can drive
             that define a survival function.",
      "i" = "{.fn ee_survival_model} takes no design matrix, so the formula
             interface cannot drive it; predict from one with
             {.fn survival_predictions}.",
      "i" = "Drop {.arg times} to predict the linear predictor instead, where
             the equation has one."
    ),
    call = NULL
  )
}

#' @noRd
abort_predict_plogit_offset <- function() {
  cli::cli_abort(
    c(
      "The offset this fit was made with cannot be carried into its predicted
       survival.",
      "i" = "{.fn plogit_predict} takes no offset, so a prediction made through
             it would silently be the prediction at an offset of zero. Writing
             the offset into the formula does not change that.",
      "i" = "Predict from the estimates directly with {.fn plogit_predict} if
             the offset is not wanted."
    ),
    call = NULL
  )
}

#' @noRd
abort_predict_plogit_exact <- function() {
  cli::cli_abort(
    c(
      "{.code deriv_method = \"exact\"} is not available for a fit of
       {.fn ee_plogit}.",
      "i" = "Its survival measure is computed through {.fn plogit_predict},
             whose matrix products cannot carry the tangents forward-mode
             differentiation works by.",
      "i" = "Use {.val capprox}, {.val fapprox}, or {.val bapprox}, or predict
             an {.fn ee_aft} fit, where {.val exact} is available."
    ),
    call = NULL
  )
}

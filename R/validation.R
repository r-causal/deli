#' Check that penalty and center have valid shapes
#'
#' Validates that `penalty` and `center` are either length 1 or the same
#' length as `theta`, and that all penalty values are non-negative.
#'
#' @param theta Numeric vector of parameters.
#' @param penalty Numeric penalty term (scalar or vector).
#' @param center Numeric center for penalty (scalar or vector).
#'
#' @returns Invisible `NULL`. Raises an error if shapes are invalid.
#' @keywords internal
check_penalty_shape <- function(theta, penalty, center) {
  # Check penalty length
  if (length(penalty) != 1 && length(penalty) != length(theta)) {
    cli::cli_abort(
      "The {.arg penalty} term must be either a single number or the same
       length as theta."
    )
  }

  # Check center length
  if (length(center) != 1 && length(center) != length(theta)) {
    cli::cli_abort(
      "The {.arg center} term must be either a single number or the same
       length as theta."
    )
  }

  # Check penalty is non-negative
  if (any(penalty < 0)) {
    cli::cli_abort("All {.arg penalty} terms must be non-negative.")
  }

  invisible(NULL)
}

#' Check that a data vector matches the observation count
#'
#' Validates that a data argument (`y`, `offset`, `q_eval`, `delta`, and the
#' like) has one value per observation, mirroring the broadcast errors Python
#' Delicatessen raises when an argument length does not match the data. Tangent
#' containers are measured on the primal, because the response reaching a
#' regression estimating equation can be a `PrimalTangentArray` whose `length()`
#' is not the observation count.
#'
#' @param x The data argument to validate.
#' @param n The number of observations the argument must match.
#' @param arg The argument name, used in the error message.
#'
#' @returns Invisible `NULL`. Raises an error if the length does not match.
#' @keywords internal
check_data_length <- function(x, n, arg) {
  len <- if (is_tangent_container(x)) {
    length(pt_arrays(x)$primal)
  } else {
    length(x)
  }
  if (len != n) {
    cli::cli_abort(
      "{.arg {arg}} must have the same length as the data ({n}), not length
       {len}."
    )
  }
  invisible(NULL)
}

#' Check that two design matrices have identical dimensions
#'
#' Validates that a counterfactual or plan design matrix has the same dimensions
#' as the observed design matrix, mirroring the explicit shape checks Python
#' Delicatessen performs in the causal estimating equations before any
#' arithmetic can silently recycle.
#'
#' @param x The reference design matrix.
#' @param y The design matrix to compare against `x`.
#' @param x_arg The reference argument name, used in the error message.
#' @param y_arg The compared argument name, used in the error message.
#'
#' @returns Invisible `NULL`. Raises an error if the dimensions differ.
#' @keywords internal
check_design_dims_match <- function(x, y, x_arg, y_arg) {
  if (!identical(dim(x), dim(y))) {
    cli::cli_abort(
      "The dimensions of {.arg {x_arg}} and {.arg {y_arg}} must be the same."
    )
  }
  invisible(NULL)
}

#' Check that a LASSO approximation epsilon is non-negative
#'
#' Validates the `epsilon` argument of the approximate LASSO estimating
#' equations directly, so a negative value is rejected with a message that names
#' `epsilon` rather than surfacing the downstream bridge-penalty message phrased
#' in terms of `gamma`. Mirrors Python Delicatessen, which validates `epsilon`
#' up front. Zero is permitted and falls through to the bridge penalty's
#' non-differentiability warning.
#'
#' @param epsilon The approximation parameter supplied by the caller.
#'
#' @returns Invisible `NULL`. Raises an error if `epsilon` is negative.
#' @keywords internal
check_epsilon <- function(epsilon) {
  if (epsilon < 0) {
    cli::cli_abort(
      "{.arg epsilon} must be greater than zero for the approximate LASSO."
    )
  }
  invisible(NULL)
}

#' Check that truncation bounds are in ascending order
#'
#' Validates that a length-2 `truncate` vector has its lower bound no
#' greater than its upper bound, mirroring the check Python Delicatessen
#' performs before clipping propensity scores.
#'
#' @param truncate Length-2 numeric vector `c(lower, upper)`.
#'
#' @returns Invisible `NULL`. Raises an error if bounds are out of order.
#' @keywords internal
check_truncate_order <- function(truncate) {
  if (truncate[1] > truncate[2]) {
    cli::cli_abort(
      "{.arg truncate} values must be specified in ascending order."
    )
  }
  invisible(NULL)
}

#' Check that survival data is valid
#'
#' Validates that event indicators are 0 or 1 (ignoring NAs) and that
#' observation times are positive (ignoring NAs).
#'
#' @param delta Numeric vector of event indicators (0 or 1).
#' @param time Numeric vector of observation times.
#'
#' @returns Invisible `NULL`. Raises an error if data is invalid.
#' @keywords internal
check_survival_data_valid <- function(delta, time) {
  # Check delta values are 0 or 1 (ignoring NAs)
  delta_no_na <- delta[!is.na(delta)]
  if (!all(delta_no_na %in% c(0, 1))) {
    cli::cli_abort(
      "All non-missing event indicator values must be either zero or one
       for survival models."
    )
  }

  # Check times are positive (ignoring NAs)
  time_no_na <- time[!is.na(time)]
  if (!all(time_no_na > 0)) {
    cli::cli_abort(
      "All non-missing observed times must be positive for survival models."
    )
  }

  invisible(NULL)
}

#' Check that init is a non-empty numeric vector
#'
#' Validates that `init` is numeric and contains at least one value.
#'
#' @param init The initial parameter vector supplied to an estimator.
#'
#' @returns Invisible `NULL`. Raises an error if `init` is invalid.
#' @keywords internal
check_estimator_init <- function(init) {
  if (!is.numeric(init)) {
    cli::cli_abort(
      "{.arg init} must be a numeric vector, not {.obj_type_friendly {init}}."
    )
  }
  if (length(init) < 1L) {
    cli::cli_abort("{.arg init} must contain at least one value.")
  }
  invisible(NULL)
}

#' Check that finite_correction names a supported correction
#'
#' Validates that `finite_correction` is either `NULL` or a supported
#' correction string. `"HC1"` is the only supported non-`NULL` value, matching
#' [finite_sample_correction()].
#'
#' @param finite_correction The finite-sample correction supplied to an
#'   estimator.
#'
#' @returns Invisible `NULL`. Raises an error if the value is unsupported.
#' @keywords internal
check_finite_correction <- function(finite_correction) {
  if (is.null(finite_correction)) {
    return(invisible(NULL))
  }
  supported <- !is.character(finite_correction) ||
    length(finite_correction) != 1L ||
    toupper(finite_correction) != "HC1"
  if (supported) {
    cli::cli_abort(
      c(
        "{.arg finite_correction} must be {.val NULL} or a supported
         correction.",
        "i" = "Supported options: {.val NULL}, {.val HC1}."
      )
    )
  }
  invisible(NULL)
}

#' Check that subset holds valid parameter indices
#'
#' Validates that `subset` is either `NULL` or a vector of whole-number,
#' 1-based parameter indices within `1:n_params`.
#'
#' @param subset The parameter subset supplied to an estimator.
#' @param n_params The number of parameters, taken from `length(init)`.
#'
#' @returns Invisible `NULL`. Raises an error if `subset` is invalid.
#' @keywords internal
check_estimator_subset <- function(subset, n_params) {
  if (is.null(subset)) {
    return(invisible(NULL))
  }
  valid <- is.numeric(subset) &&
    length(subset) >= 1L &&
    all(is.finite(subset)) &&
    all(subset == round(subset)) &&
    all(subset >= 1L) &&
    all(subset <= n_params)
  if (!valid) {
    cli::cli_abort(
      "{.arg subset} must be whole-number parameter indices between 1 and
       {n_params}."
    )
  }
  # Duplicated indices make the subsetted estimating system rank-deficient, so
  # the solver returns the starting values as if they were estimates.
  if (anyDuplicated(subset)) {
    duplicates <- unique(subset[duplicated(subset)])
    cli::cli_abort(c(
      "{.arg subset} must not contain duplicated parameter indices.",
      "i" = "Duplicated: {.val {duplicates}}."
    ))
  }
  invisible(NULL)
}

# ---- estimating-function return conditions -----------------------------------
# The aborts raised while judging what the estimating function returned carry
# condition classes, so a calling function or a test can match on the class
# rather than on the prose, as the solver warnings in `R/estimate.R` and the
# exact-mode aborts in `R/autodiff.R` do.
#
#   deli_psi_return_error
#     The estimating function returned something no sandwich can be built from.
#     Every abort in check_psi_at_init() and check_psi_at_theta() carries it, so
#     a caller who wants to recognize a refused return matches one class rather
#     than four messages. It is the family rather than the diagnosis, and the
#     narrower class below leads the vector where there is one.
#
#   deli_psi_shape_error
#     The estimating function returned a number of equations that cannot be
#     solved against the number of parameters: any mismatch under M-estimation,
#     a shortfall under GMM or at a point handed to compute_sandwich(). It is
#     the one failure raised here that a short automatic `init` can explain,
#     which is why eval_psi_at_init() catches this class alone rather than the
#     family or every error.
#
#   deli_formula_auto_init_error
#     A failure at an `init` the formula interface generated itself, reframed
#     as a failure at those automatic starting values rather than reported
#     against an argument the user never supplied. Raised only by
#     abort_formula_auto_init(), and so only on the formula path.

#' Check the estimating-function return at the initial values
#'
#' Evaluates the value of `stacked_equations(init)` and rejects returns that
#' would otherwise fail deep inside the solver with an unhelpful message: a
#' `NULL` or non-numeric return, a number of estimating equations that cannot be
#' solved against the number of parameters, or a non-finite value at the
#' starting values. This mirrors the up-front validation Python Delicatessen
#' performs before solving.
#'
#' The order the checks are judged in is load-bearing, because a single return
#' can fail more than one of them and only the first is reported. The shape is
#' two branches with different guards, an exact match under M-estimation and a
#' shortfall under GMM, and each holds its position for its own reason. Both sit
#' behind the `NULL` and numeric checks and ahead of the finite one, so a return
#' whose shape cannot be solved is reported in preference to a non-finite value
#' at the starting values, under GMM as much as under M-estimation. See the
#' comments in the body for why each check sits where it does.
#'
#' Every abort here judges the estimating function the caller supplied, from
#' several frames below the method the caller reached. This frame names the
#' parameters of an internal helper rather than anything in the call the user
#' made, so each abort reports the frame its caller passes instead.
#'
#' @param vals The value of `stacked_equations(init)`, evaluated once by the
#'   caller.
#' @param init The initial parameter vector.
#' @param allow_over_identification Logical. When `TRUE`, the number of
#'   estimating equations may exceed the number of parameters (the GMM case) and
#'   only a shortfall is rejected. Default `FALSE`.
#' @param error_call The frame to report the error against. `NULL` reports no
#'   call, which is what a caller with no frame worth naming leaves.
#'
#' @returns Invisible `NULL`. Raises an error if the return is invalid. A
#'   mismatch between the number of estimating equations and the number of
#'   parameters carries the class `deli_psi_shape_error`, which is the one
#'   failure here that an automatically generated `init` can explain. The GMM
#'   shortfall carries the number of moment conditions as the `n_moments` field
#'   of that condition as well.
#' @keywords internal
check_psi_at_init <- function(
  vals,
  init,
  allow_over_identification = FALSE,
  error_call = NULL
) {
  # `NULL` first. What it has to clear is the numeric check immediately below,
  # which is what would otherwise catch it, because `is.numeric(NULL)` is
  # `FALSE`. That check diagnoses by naming the type of the offending return,
  # and the type it names for `NULL` is `NULL`, which reports the absence of a
  # return as though it were one more type that happens not to be numeric, and
  # carries no hint; the branch here names the absence in the headline instead
  # and keeps the hint describing what a valid return holds. It has to clear
  # both shape branches too, for a separate reason: a dimensionless value
  # counts as one estimating equation there, so a `NULL` reaching them is
  # reported as a count of one against the length of `init`, which is not the
  # problem. The finite check is the only one it is unordered against, since
  # `all(is.finite(NULL))` is `TRUE`.
  if (is.null(vals)) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned {.val NULL} at the initial values.",
        "i" = "It must return a numeric vector or matrix of estimating-function
               contributions."
      ),
      class = "deli_psi_return_error",
      call = error_call
    )
  }
  # The type next, ahead of both shape branches. No length of `init`, short or
  # long, turns a numeric return into a character one, so a non-numeric return
  # is never a symptom of either count below it and counting rows on a character
  # matrix would misdirect. This also has to precede the finite check, which
  # accepts a character vector and calls every element non-finite.
  if (!is.numeric(vals)) {
    cli::cli_abort(
      "{.arg stacked_equations} must return a numeric vector or matrix at the
       initial values, not {.obj_type_friendly {vals}}.",
      class = "deli_psi_return_error",
      call = error_call
    )
  }
  # The shape ahead of the values. This is two branches under different guards,
  # never both live in one call, so their order relative to each other is free.
  # Where the pair sits relative to the other three checks is not free, and the
  # two have separate reasons for sitting here, so a maintainer moving one is
  # not free to carry the other along with it.
  n_eqs <- if (is.null(dim(vals))) 1L else nrow(vals)
  n_params <- length(init)
  # Under M-estimation, judging the values first would bury the commonest
  # length mismatch there is. An estimating function that appends a parameter
  # beyond the coefficients reads `theta[n_params + 1]`, which is `NA` when
  # `init` is one short, so its return is both the wrong number of rows and full
  # of `NA`. The row count is the accurate diagnosis of that pair, and it is the
  # one a formula fit can reframe as the automatic `init` being too short;
  # reporting the `NA`s instead names a symptom and loses the cause.
  if (!allow_over_identification && n_eqs != n_params) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned {n_eqs} estimating equation{?s} at
         the initial values, but {.arg init} has {n_params} parameter{?s}.",
        "i" = "M-estimation requires one estimating equation per parameter (a
               {n_params}-by-n matrix)."
      ),
      class = c("deli_psi_shape_error", "deli_psi_return_error"),
      call = error_call
    )
  }
  # Under GMM, a shortfall of equations is a property of the system rather than
  # of the point it was evaluated at: there is no starting value that makes an
  # under-identified system solvable. A non-finite value usually is fixable by
  # starting somewhere else, so reporting it first would send the user tuning
  # `init` against a system that has no solution to find.
  if (allow_over_identification && n_eqs < n_params) {
    cli::cli_abort(
      c(
        "The number of initial values ({n_params}) must be less than or equal to
         the number of estimating equations ({n_eqs}).",
        "i" = "GMM needs at least one moment condition per parameter, so this
               system is under-identified and has no solution at any starting
               values.",
        "i" = "Return more moment conditions from {.arg stacked_equations}, or
               estimate fewer parameters."
      ),
      class = c("deli_psi_shape_error", "deli_psi_return_error"),
      call = error_call,
      # The count travels with the condition so that the formula path can name it
      # while reframing the failure against its automatic `init`, without
      # counting the rows a second time. Only this branch carries it, so its
      # presence is also what tells that path a shortfall of moment conditions is
      # what failed rather than the M-estimation mismatch above.
      n_moments = n_eqs
    )
  }
  # Last, so it reports a non-finite value only once the return is the right
  # shape and the value is the remaining explanation.
  if (!all(is.finite(vals))) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned non-finite values at the initial
         values.",
        "i" = "This often means {.arg init} produces a divide-by-zero, a
               logarithm or square root of a non-positive number, or another
               undefined value in the estimating equations."
      ),
      class = "deli_psi_return_error",
      call = error_call
    )
  }
  invisible(NULL)
}

#' Check the estimating-function return at a point supplied as the root
#'
#' The counterpart of [check_psi_at_init()] for [compute_sandwich()], which is
#' handed a point the caller states is the root and so never runs the check
#' [estimate()] runs before solving. The returns that check rejects reach the
#' bread and the meat instead, where each either fails deep in base R against an
#' argument the caller never named or returns a matrix of the documented shape
#' built from a Jacobian that is not one.
#'
#' The checks, and the order they hold, are those of [check_psi_at_init()]; the
#' comments there give the reason each one sits where it does. Two things
#' differ. There is no `init` here and no initial values to name, so the point
#' is named by the argument that carries it. And more estimating equations than
#' parameters is the over-identified system whose sandwich this entry point is
#' asked for and returns, so only a shortfall is rejected, as it is under
#' `allow_over_identification`.
#'
#' Each abort judges an argument the caller passed to [compute_sandwich()], and
#' this function appears in no man page, so the default frame is the one frame up
#' rather than this one. `compute_sandwich()` calls this from its own body, so
#' that frame is the entry point the caller typed.
#'
#' @param vals The value of `stacked_equations(theta)`, evaluated once by the
#'   caller.
#' @param theta The parameter vector the estimating function was evaluated at.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error if the return is invalid. A
#'   shortfall of estimating equations carries the class
#'   `deli_psi_shape_error`, as the shape failures in [check_psi_at_init()] do.
#' @noRd
check_psi_at_theta <- function(vals, theta, call = rlang::caller_env()) {
  if (is.null(vals)) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned {.val NULL} at {.arg theta}.",
        "i" = "It must return a numeric vector or matrix of estimating-function
               contributions."
      ),
      class = "deli_psi_return_error",
      call = call
    )
  }
  if (!is.numeric(vals)) {
    cli::cli_abort(
      "{.arg stacked_equations} must return a numeric vector or matrix at
       {.arg theta}, not {.obj_type_friendly {vals}}.",
      class = "deli_psi_return_error",
      call = call
    )
  }
  n_eqs <- if (is.null(dim(vals))) 1L else nrow(vals)
  n_params <- length(theta)
  if (n_eqs < n_params) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned {n_eqs} estimating equation{?s} at
         {.arg theta}, which holds {n_params} parameter{?s}.",
        "i" = "The sandwich needs at least one estimating equation per
               parameter. With fewer, the bread has more columns than rows and
               its pseudo-inverse carries the shape of a covariance matrix
               without the meaning.",
        "i" = "More estimating equations than parameters is the over-identified
               system, which is accepted."
      ),
      class = c("deli_psi_shape_error", "deli_psi_return_error"),
      call = call
    )
  }
  if (!all(is.finite(vals))) {
    cli::cli_abort(
      c(
        "{.arg stacked_equations} returned non-finite values at {.arg theta}.",
        "i" = "Both the bread and the meat are built from this return, so a
               non-finite value in it leaves the whole sandwich undefined."
      ),
      class = "deli_psi_return_error",
      call = call
    )
  }
  invisible(NULL)
}

#' Evaluate and validate the estimating function at the initial values
#'
#' Performs the single evaluation of the estimating function at `init` that
#' [estimate()] needs, and validates the return with [check_psi_at_init()].
#'
#' When the formula interface generated `init` itself it records the vector on
#' the closure as the `deli_auto_init` attribute, and this reframes a failure at
#' exactly those starting values against the automatic `init`, which the user
#' never chose. Some estimating equations append parameters beyond the
#' regression coefficients (for example [ee_glm()] with `"gamma"` or
#' `"negative_binomial"`, which add a scale or dispersion parameter), so the
#' automatic `init` is one short and the estimating function either fails
#' outright or returns the wrong number of rows. Where the equation is one the
#' formula interface recognizes, the reframed message names the parameter the
#' automatic length leaves out; otherwise how much the message can claim about
#' the length depends on the failure. A wrong-shaped return is a length mismatch
#' by definition and is described as one, an error the estimating function raised
#' is described as one where the error itself reads as a length problem, and an
#' error that reads as anything else is left described as itself, so that a user
#' whose equation failed for a reason of its own is not sent looking for a length
#' that was never wrong.
#'
#' Only those two failures are reframed. A `NULL` or non-numeric return keeps
#' its own message, because an `init` one element short cannot cause either. A
#' non-finite return keeps its own message when the number of rows fits the
#' automatic length, and is reframed when it does not, because an estimating
#' function reading a parameter the automatic `init` does not reach returns
#' `NA`s and the wrong number of rows together, and the row count is the more
#' accurate of the two. Every failure raised after this point, in the solver or
#' the sandwich components, keeps its own message as well.
#'
#' @param psi The estimating-function closure.
#' @param init The initial parameter vector.
#' @param allow_over_identification Logical. When `TRUE` (the GMM case), the
#'   estimating function may return more equations than parameters, so only a
#'   shortfall is rejected. Default `FALSE`.
#' @param error_call The frame to report a failure at the caller's own `init`
#'   against, which is the [estimate()] method the caller reached. At an `init`
#'   the formula interface generated, the entry point recorded on the closure is
#'   preferred to it.
#'
#' @returns The value of `psi(init)`. Raises an error if it is not a valid
#'   estimating-function return. A failure reframed as a problem with the
#'   automatic length carries the class `deli_formula_auto_init_error`, so a
#'   caller can recognize it without matching the message.
#' @keywords internal
eval_psi_at_init <- function(
  psi,
  init,
  allow_over_identification = FALSE,
  error_call = NULL
) {
  if (!identical(attr(psi, "deli_auto_init", exact = TRUE), init)) {
    vals <- psi(init)
    check_psi_at_init(
      vals,
      init,
      allow_over_identification,
      error_call = error_call
    )
    return(vals)
  }

  # Where the estimating equation is one the formula interface recognizes, the
  # parameter the automatic length leaves out has a name, recorded beside the
  # starting values.
  appended <- attr(psi, "deli_auto_init_appended", exact = TRUE)
  # This is several frames below the method the caller reached, so the frame to
  # report the failure against travels with the starting values too. It is
  # preferred to the method for every failure at these starting values, not only
  # the ones reframed below: an `init` the caller never supplied is the formula
  # wrapper's doing, so the wrapper is the call to act on, and the whole set
  # then reports one entry point rather than dividing by which check caught the
  # failure. A closure carrying no such record falls back to the method.
  entry_point <- attr(psi, "deli_auto_init_call", exact = TRUE) %||% error_call

  # The handler covers a single call, the estimating function itself, so
  # nothing raised later can reach it.
  vals <- rlang::try_fetch(
    psi(init),
    error = function(cnd) {
      abort_formula_auto_init(
        init,
        appended,
        parent = cnd,
        error_call = entry_point
      )
    }
  )
  # Of the four returns check_psi_at_init() rejects, only a shape mismatch is a
  # length problem, so catch that class alone rather than every error. Under GMM
  # the mismatch is a shortfall of moment conditions, which is not a length
  # problem at all; the count it reports travels on the condition, so passing it
  # along is what tells the diagnostic to describe an under-identified system
  # rather than an `init` that needs lengthening.
  rlang::try_fetch(
    check_psi_at_init(
      vals,
      init,
      allow_over_identification,
      error_call = entry_point
    ),
    deli_psi_shape_error = function(cnd) {
      abort_formula_auto_init(
        init,
        appended,
        n_moments = cnd$n_moments,
        error_call = entry_point
      )
    }
  )
  vals
}

#' Abort with the automatic-`init` diagnostic
#'
#' The hint the diagnostic carries depends on what is known about the failure. A
#' GMM system that returned fewer moment conditions than there are parameters is
#' under-identified, which no `init` of any length can fix, so that failure is
#' described as itself. A recognized equation has the parameter the automatic
#' length leaves out named outright. Otherwise the length is offered as the
#' likely cause only where the failure supports it: a wrong-shaped return is a
#' length mismatch by definition, and an error the estimating function raised is
#' one when the error itself reads as a length problem. An error that reads as
#' anything else gets a hint that describes only what is known, because a user
#' whose equation failed for a reason of its own was sent looking for a length
#' that was never wrong.
#'
#' @param init The automatically generated initial parameter vector.
#' @param appended The parameter the estimating equation estimates beyond the
#'   design coefficients, or `NULL` when the equation is not one this package
#'   recognizes.
#' @param n_moments The number of moment conditions a GMM system returned, when
#'   that count fell short of the number of parameters, and `NULL` otherwise.
#' @param parent The error raised while evaluating the estimating function, or
#'   `NULL` when the estimating function returned a wrong-shaped value instead
#'   of failing.
#' @param error_call The frame to report the error against, recorded on the
#'   estimating-function closure by `prepare_formula_psi()`. `NULL` reports no
#'   call, which is what a closure carrying no such record leaves.
#'
#' @returns Never returns; always raises an error carrying the class
#'   `deli_formula_auto_init_error`.
#' @noRd
abort_formula_auto_init <- function(
  init,
  appended = NULL,
  n_moments = NULL,
  parent = NULL,
  error_call = NULL
) {
  n_params <- length(init)
  # An error means the estimating function did not even evaluate at the
  # automatic init, so its cause is unknown; do not assert it is a length
  # problem. A clean return that simply has the wrong shape is a genuine length
  # mismatch and is described as such.
  errored <- !is.null(parent)
  header <- if (errored) {
    "Evaluating the estimating function at the automatic zero {.arg init} of
     length {n_params} (the number of model-matrix columns) failed."
  } else {
    "The automatic zero {.arg init} has length {n_params}, the number of
     model-matrix columns, which does not fit the estimating function."
  }
  length_hint <- "Estimating equations such as {.fn ee_glm} with {.val gamma} or
     {.val negative_binomial} append an extra parameter and need an
     {.arg init} one longer than the coefficients."
  bullets <- if (!is.null(n_moments)) {
    # An under-identified GMM system is judged first, ahead of everything the
    # automatic length could explain, because it is a property of the system: no
    # `init` of any length has a solution to find, and a longer one only widens
    # the shortfall. The two counts stay, since they are what the reader checks
    # the diagnosis against.
    c(
      "i" = "The estimating function returned {n_moments} moment
             condition{?s} for {n_params} parameter{?s}, so the system is
             under-identified. GMM needs at least one moment condition per
             parameter.",
      "i" = "Return more moment conditions from {.arg .ee}, or fit a formula
             with fewer terms."
    )
  } else if (!is.null(appended)) {
    # A recognized equation can have the missing parameter named and counted
    # rather than described, on both paths, and whatever the failure was. It is
    # the equation's own layout that is known here, not anything read off the
    # error, so this branch is judged first.
    c(
      "i" = "This estimating equation estimates one parameter beyond the
             design coefficients, {.val {appended}}, so it needs an
             {.arg init} of length {n_params + 1}.",
      "i" = "Supply an explicit {.arg init} of length {n_params + 1}."
    )
  } else if (errored && !error_looks_length_related(parent)) {
    # Nothing is known about an unrecognized equation's parameters, and an
    # error that does not read as a length problem is no evidence that the
    # length is one. Both bullets stay with what is known: the equation failed,
    # its own error says why, and an explicit `init` is how to rule the
    # automatic one out rather than the fix being asserted in advance.
    c(
      "i" = "The estimating function itself failed at those starting values;
             the error it raised is reported below.",
      "i" = "Supplying an explicit {.arg init} rules the automatic one out as
             the cause."
    )
  } else {
    # A wrong-shaped return, or an error that reads as a length problem. Only
    # the second of those needs the cause offered as likely rather than stated,
    # since the first is a length mismatch on the face of it.
    hint <- if (errored) {
      paste("A length mismatch is the most common cause.", length_hint)
    } else {
      length_hint
    }
    c("i" = hint, "i" = "Supply an explicit {.arg init} of the correct length.")
  }
  cli::cli_abort(
    c(header, bullets),
    parent = parent,
    call = error_call,
    class = "deli_formula_auto_init_error"
  )
}

#' Whether an error reads as a length problem
#'
#' Decides whether the error an estimating function raised at the automatic
#' `init` supports naming the length as the likely cause. The families a
#' too-short automatic `init` actually produces are three: a matrix product
#' against a parameter vector one element longer than the parameters reach,
#' which R reports as non-conformable; `[[` past the end of the parameter
#' vector; and the equation's own check on the length of `theta`, which says so
#' in its own words.
#'
#' Out-of-bounds `[[` carries a condition class of its own, so that family is
#' recognized without reading a message at all. The other two are matched on the
#' text, which two things follow from. The first is that the parameter vector's
#' name is not part of the match: the offending call belongs to whatever the
#' equation calls its parameters, and [ee_glm()] with `"gamma"` fails in
#' `X %*% beta`, so a match requiring the name `theta` would miss the very case
#' this diagnostic was written for. The second is that a translated message will
#' not match, which costs the hint rather than misdirecting: the caller falls
#' back to the neutral wording, and the error itself is reported either way.
#'
#' @param cnd The condition raised while evaluating the estimating function.
#'
#' @returns `TRUE` when the error reads as a length problem, `FALSE` otherwise.
#' @noRd
error_looks_length_related <- function(cnd) {
  if (inherits(cnd, "subscriptOutOfBoundsError")) {
    return(TRUE)
  }
  any(grepl(
    "non-conformable|subscript out of bounds|length",
    conditionMessage(cnd)
  ))
}

# ---- formula-interface conditions --------------------------------------------
# The conditions the formula interface raises about the estimating equation it
# was handed carry a condition class, so a caller or a test can match on the
# class rather than on the prose. Each of them reports the entry point the caller
# typed, `m_estimate()` or `gmm_estimate()`, which the `.formula` methods pass
# down because every frame these are raised in belongs to a helper no caller
# wrote.
#
#   deli_formula_ee_signature_error
#     The arguments of `.ee` leave something the interface fills nowhere to go:
#     `theta`, the model matrix it passes as `X`, the response it passes
#     positionally, or the offset an `offset()` term in the formula supplies.
#     Raised by check_formula_ee_signature().
#
#   deli_formula_ee_argument_error
#     A name the caller supplied in `...` matches no argument of `.ee` exactly.
#     Raised by check_formula_ee_dots().
#
#   deli_formula_auto_init_error
#     The automatic `init` does not fit the estimating equation. Raised by
#     abort_formula_auto_init() above, at the point the estimating function is
#     first evaluated rather than where the formula was read.

#' Check that `.ee` can receive what the formula interface fills
#'
#' The formula interface fills four arguments of its own: `theta`, the model
#' matrix it passes as `X`, the response it passes positionally, and the offset
#' it takes from an `offset()` term in the formula. An equation whose arguments
#' leave any of them nowhere to go cannot be driven by a formula at all, and
#' reaching it anyway produced base R's unused-argument error, which pastes the
#' whole offending vector into its message, so a design matrix or a response was
#' reported one line per observation long.
#'
#' Runs before the exact match on the names in `...`, because a name the caller
#' supplied is worth reporting only once the equation can be driven at all. The
#' two checks divide by whose argument is at fault: the arguments the interface
#' fills here, the names the caller wrote there.
#'
#' @param .ee The estimating-equation function passed to the formula interface.
#' @param ee_args The evaluated `...` arguments forwarded to it, after any
#'   `offset()` term in the formula has been added to them.
#' @param has_formula_offset Logical. Did the formula carry an `offset()` term?
#' @param ee_name The name `.ee` was passed under, or `NULL` when it arrived as
#'   an anonymous function.
#' @param error_call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_formula_ee_signature_error` if an argument the interface fills has
#'   nowhere to go.
#' @noRd
check_formula_ee_signature <- function(
  .ee,
  ee_args,
  has_formula_offset = FALSE,
  ee_name = NULL,
  error_call = NULL
) {
  formal_names <- formula_ee_formals(.ee)
  if (is.null(formal_names)) {
    return(invisible(NULL))
  }
  # An equation with a `...` of its own has somewhere to put every argument,
  # so nothing the interface fills can be left out.
  if ("..." %in% formal_names) {
    return(invisible(NULL))
  }

  supplied <- formula_ee_dots_names(ee_args)

  faults <- character()
  if (!"theta" %in% formal_names) {
    faults <- c(
      faults,
      "x" = "It has no {.arg theta} argument for the parameter vector."
    )
  }
  if (!"X" %in% formal_names) {
    faults <- c(
      faults,
      "x" = "It has no {.arg X} argument for the model matrix."
    )
  }
  if (is.null(formula_ee_response_slot(formal_names, supplied))) {
    faults <- c(
      faults,
      "x" = "It has no argument left for the response, which is passed
             positionally after {.arg theta} and {.arg X}."
    )
  }

  if (length(faults) > 0L) {
    equation <- formula_ee_label(ee_name)
    cli::cli_abort(
      c(
        "{equation} cannot be driven by a formula.",
        faults,
        # An equation with no arguments at all has nothing to list.
        if (length(formal_names) > 0L) {
          c("i" = "It takes {.arg {formal_names}}.")
        },
        "i" = "Fit it through the function interface instead, by passing a
               function of {.arg theta} as {.arg stacked_equations}."
      ),
      call = error_call,
      class = "deli_formula_ee_signature_error"
    )
  }

  # The offset is the one argument the interface fills from the formula rather
  # than from the model frame, so the report names the term it came from. The
  # caller wrote no `offset` name at all, so nothing here is a misspelling that a
  # correction could be suggested for.
  if (has_formula_offset && !"offset" %in% formal_names) {
    equation <- formula_ee_label(ee_name)
    cli::cli_abort(
      c(
        "The {.code offset()} term in the formula has no argument to go to.",
        "x" = "{equation} has no {.arg offset} argument, and an {.code offset()}
               term reaches the estimating equation through one.",
        "i" = "Drop the term from the formula, or fit an equation that takes an
               offset."
      ),
      call = error_call,
      class = "deli_formula_ee_signature_error"
    )
  }

  invisible(NULL)
}

#' Match the names supplied in `...` against the arguments of `.ee` exactly
#'
#' R matches a supplied name that is a prefix of exactly one formal to that
#' formal, and no built-in estimating equation takes a `...` of its own for
#' another name to fall into, so forwarding `...` with [do.call()] resolved an
#' abbreviation or a prefix typo against the equation's arguments: `weight = w`
#' reached `ee_regression()`'s `weights` and returned the weighted estimates.
#' Such a fit succeeded silently while reporting different numbers, so the names
#' are matched exactly here and anything else is refused before the equation is
#' evaluated.
#'
#' @param .ee The estimating-equation function passed to the formula interface.
#' @param ee_args The evaluated `...` arguments forwarded to it.
#' @param ee_name The name `.ee` was passed under, or `NULL` when it arrived as
#'   an anonymous function.
#' @param error_call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_formula_ee_argument_error` if a supplied name matches no argument.
#' @noRd
check_formula_ee_dots <- function(
  .ee,
  ee_args,
  ee_name = NULL,
  error_call = NULL
) {
  formal_names <- formula_ee_formals(.ee)
  if (is.null(formal_names)) {
    return(invisible(NULL))
  }
  # A `...` of the equation's own absorbs any name, so the caller may write one
  # that no argument of it accounts for.
  if ("..." %in% formal_names) {
    return(invisible(NULL))
  }

  supplied <- formula_ee_dots_names(ee_args)
  refused <- setdiff(supplied, formal_names)
  if (length(refused) == 0L) {
    return(invisible(NULL))
  }

  equation <- formula_ee_label(ee_name)

  # The arguments a caller may write are the ones left once the interface has
  # filled its own, so neither the list nor the suggestion sends them after an
  # argument they are not allowed to pass.
  reserved <- c(
    "theta",
    "X",
    formula_ee_response_slot(formal_names, supplied)
  )
  candidates <- setdiff(formal_names, reserved)
  # One refused name can be answered with the argument it was meant for. Several
  # cannot: a suggestion for one of them reads as a suggestion for all.
  meant <- if (length(refused) == 1L) {
    formula_ee_suggestion(refused, candidates)
  } else {
    character()
  }

  cli::cli_abort(
    c(
      "{equation} has no {.arg {refused}} argument{?s}.",
      if (length(meant) > 0L) {
        c("i" = "Did you mean {.or {.arg {meant}}}?")
      } else if (length(candidates) > 0L) {
        c(
          "i" = "Besides the arguments the formula interface fills itself, it
                 takes {.arg {candidates}}."
        )
      } else {
        c(
          "i" = "It takes no argument the formula interface does not fill itself."
        )
      }
    ),
    call = error_call,
    class = "deli_formula_ee_argument_error"
  )
}

#' The argument names of an estimating equation, where they can be read
#'
#' `args()` gives an argument list for anything that has one and `NULL` for
#' everything else: a primitive such as `[` has no argument list of its own, and
#' neither has a value that is not a function at all. `formals()` warns rather
#' than reporting anything on those, so `NULL` here means the arguments cannot be
#' read, and the checks decline to speak about a function whose arguments they
#' cannot read. Such an `.ee` fails when it is called instead. An equation that
#' takes no arguments is a readable, empty list rather than an unreadable one.
#'
#' @param .ee The estimating-equation function passed to the formula interface.
#'
#' @returns A character vector of argument names, or `NULL` when they cannot be
#'   read.
#' @noRd
formula_ee_formals <- function(.ee) {
  fn <- args(.ee)
  if (!is.function(fn)) {
    return(NULL)
  }
  names(formals(fn)) %||% character()
}

#' The argument the formula response is passed to
#'
#' The response is passed positionally, so it fills the first argument that is
#' neither `theta` nor `X` nor matched by a name the caller supplied. `NULL` when
#' the equation has no such argument, which is what leaves the response nowhere
#' to go.
#'
#' @param formal_names The argument names of the estimating equation.
#' @param supplied The names the caller supplied in `...`.
#'
#' @returns A string, or `NULL`.
#' @noRd
formula_ee_response_slot <- function(formal_names, supplied) {
  free <- setdiff(formal_names, c("theta", "X", supplied))
  if (length(free) == 0L) NULL else free[[1]]
}

#' The names a caller supplied in `...`
#'
#' Drops the empty name an unnamed argument carries, which is passed positionally
#' and so matches no argument by name at all.
#'
#' @param ee_args The evaluated `...` arguments.
#'
#' @returns A character vector, empty when nothing was named.
#' @noRd
formula_ee_dots_names <- function(ee_args) {
  supplied <- names(ee_args)
  if (is.null(supplied)) {
    return(character())
  }
  supplied[nzchar(supplied)]
}

#' The argument a refused name was probably meant to be
#'
#' A name that prefixes an argument is answered with that argument exactly,
#' because R would have matched it there and the fit it produced is the one the
#' caller is trying to explain. A name that prefixes nothing is answered from its
#' spelling, within two edits, which covers a transposition and a doubled or
#' dropped pair of letters. Anything further away is left unanswered: a
#' suggestion nothing supports sends the caller after the wrong argument.
#'
#' @param name The refused name.
#' @param candidates The arguments the caller may pass.
#'
#' @returns A character vector of candidates, empty when none is close enough.
#' @noRd
formula_ee_suggestion <- function(name, candidates) {
  if (length(candidates) == 0L) {
    return(character())
  }
  prefixed <- candidates[startsWith(candidates, name)]
  if (length(prefixed) > 0L) {
    return(prefixed)
  }
  distances <- as.integer(utils::adist(name, candidates, ignore.case = TRUE))
  candidates[distances <= 2L]
}

#' How to name the estimating equation in a report about it
#'
#' An equation passed as a name is named, which says which of the arguments in
#' the call is the one to change. An anonymous function has no name to report and
#' is described by the argument it arrived in.
#'
#' @param ee_name The name `.ee` was passed under, or `NULL`.
#'
#' @returns A formatted string for interpolation into a cli message.
#' @noRd
formula_ee_label <- function(ee_name) {
  if (is.null(ee_name)) {
    cli::format_inline("The estimating equation passed as {.arg .ee}")
  } else {
    cli::format_inline("{.fn {ee_name}}")
  }
}

#' Check the return of a custom solver
#'
#' Validates that a user-supplied solver returned a numeric vector of the
#' expected length, so a malformed return produces an informative error rather
#' than an opaque failure while assembling the sandwich components.
#'
#' @param theta The value returned by the custom solver.
#' @param n_params The expected number of solved parameters.
#'
#' @returns Invisible `NULL`. Raises an error if the return is invalid.
#' @keywords internal
check_solver_return <- function(theta, n_params) {
  if (is.null(theta) || !is.numeric(theta) || length(theta) != n_params) {
    cli::cli_abort(c(
      "The custom {.arg solver} must return a numeric vector of length
       {n_params}.",
      "i" = "It returned {.obj_type_friendly {theta}} of length
             {length(theta)}."
    ))
  }
  invisible(NULL)
}

#' Normalize and validate a derivative-method argument
#'
#' Lowercases `deriv_method` and validates it against the supported options,
#' mirroring Python Delicatessen, which lowercases every method comparison and
#' accepts any case. Returns the normalized value so callers can branch on it
#' directly.
#'
#' @param deriv_method The derivative method supplied by the caller. One of
#'   `"capprox"`, `"fapprox"`, `"bapprox"`, or `"exact"`, in any case.
#' @param call The frame to report the refusal against. This function judges an
#'   argument the caller wrote and appears in no man page, so the default is the
#'   frame one up rather than this one. Every caller runs the check in its own
#'   body, so that frame is the entry point the caller reached, except in
#'   `delta_method_impl()`, which is itself a worker and passes the frame it was
#'   given.
#'
#' @returns The lower-cased method string. Raises an error if the value is not a
#'   single supported string.
#' @keywords internal
check_deriv_method <- function(deriv_method, call = rlang::caller_env()) {
  if (!is.character(deriv_method) || length(deriv_method) != 1L) {
    cli::cli_abort(
      "{.arg deriv_method} must be a single string, not
       {.obj_type_friendly {deriv_method}}.",
      call = call
    )
  }
  normalized <- tolower(deriv_method)
  supported <- c("capprox", "fapprox", "bapprox", "exact")
  if (!normalized %in% supported) {
    cli::cli_abort(
      c(
        "{.arg deriv_method} value {.val {deriv_method}} is not supported.",
        "i" = "Supported options: {.val capprox}, {.val fapprox},
               {.val bapprox}, {.val exact}."
      ),
      call = call
    )
  }
  normalized
}

#' Check that a distribution or link is a single string
#'
#' The generalized linear estimating equations read `distribution` and `link`
#' by name, and the helpers that read them dispatch through `if` on a comparison
#' with a single value. A `NULL` or a vector longer than one never reaches the
#' unsupported-name diagnostic those helpers raise: it fails the `if` first, as
#' base R's `argument is of length zero` or `the condition has length > 1`,
#' reported against a branch the caller never wrote. `ee_glm()` partitions
#' `theta` on a comparison of its own before either helper is reached, so the
#' name is judged where the caller supplied it instead.
#'
#' @param value The distribution or link supplied by the caller.
#' @param arg The argument name, used in the error message.
#' @param call The frame to report the refusal against, which is the estimating
#'   equation the caller wrote.
#'
#' @returns Invisible `NULL`. Raises an error if the value is not a single
#'   non-missing string.
#' @keywords internal
check_family_name <- function(value, arg, call = rlang::caller_env()) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    "{.arg {arg}} must be a single string, not {.obj_type_friendly {value}}.",
    call = call
  )
}

#' Check that a finite-difference step is a single positive number
#'
#' Validates `dx`, the absolute perturbation the finite-difference methods
#' apply. Nothing downstream reports a step that cannot be taken, because
#' [approx_differentiation()] floors each parameter's step at that magnitude's
#' floating-point resolution and the floor absorbs the two values a caller is
#' most likely to mean something by. A `dx` of zero or a negative `dx` becomes
#' the floor at every parameter away from zero, silently substituting a step the
#' caller did not ask for, and leaves a division by zero at a parameter of
#' exactly zero, where the floor is zero as well. A longer vector is recycled one
#' element per parameter, so each parameter is differentiated with a different
#' step and no two rows of the Jacobian are comparable.
#'
#' The step is validated wherever it is supplied, including under
#' `deriv_method = "exact"` and on a prediction that asks for no standard error,
#' neither of which takes a step at all. A value that cannot be a step is worth
#' reporting whether or not this particular call would have used it, since the
#' alternative is accepting it in silence and rejecting it on the next call.
#'
#' @param dx The finite-difference step supplied by the caller.
#' @param call The frame to report the refusal against, on the same terms as
#'   [check_deriv_method()]'s.
#'
#' @returns Invisible `NULL`. Raises an error if the step is not a single
#'   positive finite number.
#' @keywords internal
check_dx <- function(dx, call = rlang::caller_env()) {
  if (is.numeric(dx) && length(dx) == 1L && is.finite(dx) && dx > 0) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "{.arg dx} must be a single positive finite number.",
      # A single number is reported by value, because naming its type reports it
      # as "a number", which says nothing about which of the three requirements
      # it missed. Anything else is reported by type, which is what it failed on.
      "x" = if (is.numeric(dx) && length(dx) == 1L) {
        "The step supplied was {.val {dx}}."
      } else {
        "The step supplied was {.obj_type_friendly {dx}}."
      }
    ),
    call = call
  )
}

#' Check that an over-identification control is a single number
#'
#' Validates that a GMM over-identification control (`overid_maxiter` or
#' `overid_tolerance`) is a single, non-missing number. Non-positive values are
#' permitted: they encode degenerate but supported settings that match Python
#' Delicatessen, so they are not rejected here.
#'
#' @param value The control value supplied to [GMMEstimator()].
#' @param arg The argument name, used in the error message.
#'
#' @returns Invisible `NULL`. Raises an error if the value is invalid.
#' @keywords internal
check_overid_scalar <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    cli::cli_abort("{.arg {arg}} must be a single number.")
  }
  invisible(NULL)
}

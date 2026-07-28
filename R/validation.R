#' Check that penalty and center have valid shapes
#'
#' Validates that `penalty` and `center` are either length 1 or the same
#' length as `theta`, and that all penalty values are non-negative.
#'
#' @param theta Numeric vector of parameters.
#' @param penalty Numeric penalty term (scalar or vector).
#' @param center Numeric center for penalty (scalar or vector).
#'
#' @return Invisible `NULL`. Raises an error if shapes are invalid.
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
#' @return Invisible `NULL`. Raises an error if the length does not match.
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
#' @return Invisible `NULL`. Raises an error if the dimensions differ.
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
#' @return Invisible `NULL`. Raises an error if `epsilon` is negative.
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
#' @return Invisible `NULL`. Raises an error if bounds are out of order.
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
#' @return Invisible `NULL`. Raises an error if data is invalid.
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
#' @return Invisible `NULL`. Raises an error if `init` is invalid.
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
#' @return Invisible `NULL`. Raises an error if the value is unsupported.
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
#' @return Invisible `NULL`. Raises an error if `subset` is invalid.
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

# ---- initial-value failure conditions ---------------------------------------
# The aborts raised while judging the estimating function at the initial values
# carry condition classes, so a calling function or a test can match on the
# class rather than on the prose, as the solver warnings in `R/estimate.R` and
# the exact-mode aborts in `R/autodiff.R` do.
#
#   deli_psi_shape_error
#     The estimating function returned a number of equations that cannot be
#     solved against the number of parameters: any mismatch under M-estimation,
#     a shortfall under GMM. It is the one failure raised here that a short
#     automatic `init` can explain, which is why eval_psi_at_init() catches
#     this class alone rather than every error.
#
#   deli_formula_auto_init_error
#     A failure at an `init` the formula interface generated itself, reframed
#     as a problem with the automatic length rather than reported against an
#     argument the user never supplied. Raised only by
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
#' @param vals The value of `stacked_equations(init)`, evaluated once by the
#'   caller.
#' @param init The initial parameter vector.
#' @param allow_over_identification Logical. When `TRUE`, the number of
#'   estimating equations may exceed the number of parameters (the GMM case) and
#'   only a shortfall is rejected. Default `FALSE`.
#'
#' @return Invisible `NULL`. Raises an error if the return is invalid. A
#'   mismatch between the number of estimating equations and the number of
#'   parameters carries the class `deli_psi_shape_error`, which is the one
#'   failure here that an automatically generated `init` can explain.
#' @keywords internal
check_psi_at_init <- function(vals, init, allow_over_identification = FALSE) {
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
    cli::cli_abort(c(
      "{.arg stacked_equations} returned {.val NULL} at the initial values.",
      "i" = "It must return a numeric vector or matrix of estimating-function
             contributions."
    ))
  }
  # The type next, ahead of both shape branches. No length of `init`, short or
  # long, turns a numeric return into a character one, so a non-numeric return
  # is never a symptom of either count below it and counting rows on a character
  # matrix would misdirect. This also has to precede the finite check, which
  # accepts a character vector and calls every element non-finite.
  if (!is.numeric(vals)) {
    cli::cli_abort(
      "{.arg stacked_equations} must return a numeric vector or matrix at the
       initial values, not {.obj_type_friendly {vals}}."
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
      class = "deli_psi_shape_error"
    )
  }
  # Under GMM, a shortfall of equations is a property of the system rather than
  # of the point it was evaluated at: there is no starting value that makes an
  # under-identified system solvable. A non-finite value usually is fixable by
  # starting somewhere else, so reporting it first would send the user tuning
  # `init` against a system that has no solution to find.
  if (allow_over_identification && n_eqs < n_params) {
    cli::cli_abort(
      "The number of initial values ({n_params}) must be less than or equal to
       the number of estimating equations ({n_eqs}).",
      class = "deli_psi_shape_error"
    )
  }
  # Last, so it reports a non-finite value only once the return is the right
  # shape and the value is the remaining explanation.
  if (!all(is.finite(vals))) {
    cli::cli_abort(c(
      "{.arg stacked_equations} returned non-finite values at the initial
       values.",
      "i" = "This often means {.arg init} produces a divide-by-zero, a
             logarithm or square root of a non-positive number, or another
             undefined value in the estimating equations."
    ))
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
#' exactly those starting values as a problem with the automatic length, which
#' the user never chose. Some estimating equations append parameters beyond the
#' regression coefficients (for example [ee_glm()] with `"gamma"` or
#' `"negative_binomial"`, which add a scale or dispersion parameter), so the
#' automatic `init` is one short and the estimating function either fails
#' outright or returns the wrong number of rows.
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
#'
#' @return The value of `psi(init)`. Raises an error if it is not a valid
#'   estimating-function return. A failure reframed as a problem with the
#'   automatic length carries the class `deli_formula_auto_init_error`, so a
#'   caller can recognize it without matching the message.
#' @keywords internal
eval_psi_at_init <- function(psi, init, allow_over_identification = FALSE) {
  if (!identical(attr(psi, "deli_auto_init", exact = TRUE), init)) {
    vals <- psi(init)
    check_psi_at_init(vals, init, allow_over_identification)
    return(vals)
  }

  # The handler covers a single call, the estimating function itself, so
  # nothing raised later can reach it.
  vals <- rlang::try_fetch(
    psi(init),
    error = function(cnd) abort_formula_auto_init(init, parent = cnd)
  )
  # Of the four returns check_psi_at_init() rejects, only a shape mismatch is a
  # length problem, so catch that class alone rather than every error.
  rlang::try_fetch(
    check_psi_at_init(vals, init, allow_over_identification),
    deli_psi_shape_error = function(cnd) abort_formula_auto_init(init)
  )
  vals
}

#' Abort with the automatic-`init` diagnostic
#'
#' @param init The automatically generated initial parameter vector.
#' @param parent The error raised while evaluating the estimating function, or
#'   `NULL` when the estimating function returned a wrong-shaped value instead
#'   of failing.
#'
#' @return Never returns; always raises an error carrying the class
#'   `deli_formula_auto_init_error`.
#' @noRd
abort_formula_auto_init <- function(init, parent = NULL) {
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
  hint <- if (errored) {
    "A length mismatch is the most common cause. Estimating equations such as
     {.fn ee_glm} with {.val gamma} or {.val negative_binomial} append an extra
     parameter and need an {.arg init} one longer than the coefficients."
  } else {
    "Estimating equations such as {.fn ee_glm} with {.val gamma} or
     {.val negative_binomial} append an extra parameter and need an
     {.arg init} one longer than the coefficients."
  }
  cli::cli_abort(
    c(
      header,
      "i" = hint,
      "i" = "Supply an explicit {.arg init} of the correct length."
    ),
    parent = parent,
    call = NULL,
    class = "deli_formula_auto_init_error"
  )
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
#' @return Invisible `NULL`. Raises an error if the return is invalid.
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
#'
#' @return The lower-cased method string. Raises an error if the value is not a
#'   single supported string.
#' @keywords internal
check_deriv_method <- function(deriv_method) {
  if (!is.character(deriv_method) || length(deriv_method) != 1L) {
    cli::cli_abort(
      "{.arg deriv_method} must be a single string, not
       {.obj_type_friendly {deriv_method}}."
    )
  }
  normalized <- tolower(deriv_method)
  supported <- c("capprox", "fapprox", "bapprox", "exact")
  if (!normalized %in% supported) {
    cli::cli_abort(c(
      "{.arg deriv_method} value {.val {deriv_method}} is not supported.",
      "i" = "Supported options: {.val capprox}, {.val fapprox},
             {.val bapprox}, {.val exact}."
    ))
  }
  normalized
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
#' @return Invisible `NULL`. Raises an error if the value is invalid.
#' @keywords internal
check_overid_scalar <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    cli::cli_abort("{.arg {arg}} must be a single number.")
  }
  invisible(NULL)
}

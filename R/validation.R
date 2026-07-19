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
      "All non-missing observed times must be non-negative for survival
       models."
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
      "{.arg init} must be a numeric vector, not {.obj_type_of {init}}."
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

#' Check the estimating-function return at the initial values
#'
#' Evaluates the value of `stacked_equations(init)` and rejects returns that
#' would otherwise fail deep inside the solver with an unhelpful message: a
#' `NULL` or non-numeric return, a non-finite value at the starting values, or
#' (for M-estimation) a number of estimating equations that does not match the
#' number of parameters. This mirrors the up-front validation Python
#' Delicatessen performs before solving.
#'
#' @param vals The value of `stacked_equations(init)`, evaluated once by the
#'   caller.
#' @param init The initial parameter vector.
#' @param allow_over_identification Logical. When `TRUE`, the number of
#'   estimating equations may exceed the number of parameters (the GMM case) and
#'   the dimension check is skipped. Default `FALSE`.
#'
#' @return Invisible `NULL`. Raises an error if the return is invalid.
#' @keywords internal
check_psi_at_init <- function(vals, init, allow_over_identification = FALSE) {
  if (is.null(vals)) {
    cli::cli_abort(c(
      "{.arg stacked_equations} returned {.val NULL} at the initial values.",
      "i" = "It must return a numeric vector or matrix of estimating-function
             contributions."
    ))
  }
  if (!is.numeric(vals)) {
    cli::cli_abort(
      "{.arg stacked_equations} must return a numeric vector or matrix at the
       initial values, not {.obj_type_of {vals}}."
    )
  }
  if (!all(is.finite(vals))) {
    cli::cli_abort(c(
      "{.arg stacked_equations} returned non-finite values at the initial
       values.",
      "i" = "This often means {.arg init} produces a divide-by-zero, a
             logarithm or square root of a non-positive number, or another
             undefined value in the estimating equations."
    ))
  }
  if (!allow_over_identification) {
    n_eqs <- if (is.null(dim(vals))) 1L else nrow(vals)
    n_params <- length(init)
    if (n_eqs != n_params) {
      cli::cli_abort(c(
        "{.arg stacked_equations} returned {n_eqs} estimating equation{?s} at
         the initial values, but {.arg init} has {n_params} parameter{?s}.",
        "i" = "M-estimation requires one estimating equation per parameter (a
               {n_params}-by-n matrix)."
      ))
    }
  }
  invisible(NULL)
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
      "i" = "It returned {.obj_type_of {theta}} of length {length(theta)}."
    ))
  }
  invisible(NULL)
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

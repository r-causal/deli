#' Generate observation weights
#'
#' Returns observation weights. If `weights` is `NULL`, returns a vector of ones
#' (equal weighting). A numeric vector is validated against `n` and returned. A
#' weight matrix (one row per observation, one column per time interval, used by
#' [ee_plogit()] for time-varying weights) is validated on its row count and
#' returned with its dimensions preserved; the column-count check is left to the
#' caller, which alone knows the number of intervals.
#'
#' @param n Integer number of observations.
#' @param weights Numeric vector of weights, an n-row weight matrix, or `NULL`
#'   for equal weights.
#'
#' @returns A numeric vector of length `n`, or the supplied n-row matrix.
#'
#' @keywords internal
generate_weights <- function(n, weights = NULL) {
  if (is.null(weights)) {
    return(rep(1, n))
  }

  # A weight matrix carries one row per observation. Validate the row count and
  # return the matrix untouched; as.numeric() would strip the dimensions that
  # the caller needs to broadcast weights across time intervals.
  if (is.matrix(weights)) {
    if (nrow(weights) != n) {
      cli::cli_abort(
        "{.arg weights} must have one row per observation ({n}),
         not {nrow(weights)} rows."
      )
    }
    return(weights)
  }

  # A tangent-carrying weights vector arises under exact-mode autodiff: the
  # weights are derived from theta, as in ee_ipw_msm(), which forms IPW weights
  # from the propensity model and passes them to a weighted GLM. as.numeric()
  # would strip the primal-tangent structure and abort, so validate the length on
  # the primal and pass the tangent-carrying vector through (normalized so a
  # downstream w * score differentiates). This mirrors the is.matrix
  # short-circuit above; the numeric path below is left bit-identical.
  if (is_tangent_container(weights)) {
    primal_length <- length(pt_arrays(weights)$primal)
    if (primal_length != n) {
      cli::cli_abort(
        "{.arg weights} must have the same length as the data ({n}),
         not length {primal_length}."
      )
    }
    return(pt_as_vector(weights))
  }

  if (length(weights) != n) {
    cli::cli_abort(
      "{.arg weights} must have the same length as the data ({n}),
       not length {length(weights)}."
    )
  }

  as.numeric(weights)
}

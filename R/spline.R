#' Generate polynomial spline basis terms
#'
#' Generates polynomial spline terms for a numeric vector at pre-specified
#' knot locations. Default is restricted (natural) cubic splines, but
#' unrestricted splines with different polynomial terms can also be generated.
#'
#' Unrestricted splines for knot \eqn{k} are:
#' \deqn{s_k(X) = I(X > k) (X - k)^a}
#'
#' Restricted (natural) splines subtract the last spline term:
#' \deqn{r_k(X) = s_k(X) - s_K(X)}
#' where \eqn{K} is the largest knot. Restricted splines return one fewer
#' column than the number of knots.
#'
#' @param x Numeric vector of observed values.
#' @param knots Numeric vector of knot locations. Should be between the
#'   min and max of `x`.
#' @param power Numeric power for the spline terms. Default `3` (cubic).
#' @param restricted Logical. If `TRUE` (default), generate restricted
#'   (natural) splines. If `FALSE`, generate unrestricted splines.
#' @param normalized Logical. If `TRUE`, divide spline terms by the range
#'   of knots raised to `power`. Default `FALSE`.
#'
#' @returns A matrix with `length(x)` rows. Number of columns is
#'   `length(knots)` for unrestricted or `length(knots) - 1` for restricted.
#'
#' @export
deli_spline <- function(
  x,
  knots,
  power = 3,
  restricted = TRUE,
  normalized = FALSE
) {
  # Sort knots and set up output

  knots <- sort(as.numeric(unlist(knots)))
  n_knots <- length(knots)
  x <- as.numeric(x)
  n_obs <- length(x)

  # Create the spline basis matrix

  spline_terms <- matrix(0, nrow = n_obs, ncol = n_knots)

  # Determine normalization divisor

  if (normalized) {
    if (n_knots == 1L) {
      divisor <- knots[1]^power
    } else {
      divisor <- (knots[n_knots] - knots[1])^power
    }
  } else {
    divisor <- 1
  }

  # Generate each spline term

  for (i in seq_len(n_knots)) {
    # Truncated power basis: (x - knot)^power when x > knot, else 0
    vals <- ifelse(x > knots[i], (x - knots[i])^power, 0)
    # Preserve NAs
    spline_terms[, i] <- ifelse(is.na(x), NA_real_, vals)
  }

  # Apply restriction (natural spline) if requested

  if (restricted) {
    for (i in seq_len(n_knots - 1L)) {
      vals <- ifelse(
        x > knots[i],
        spline_terms[, i] - spline_terms[, n_knots],
        0
      )
      spline_terms[, i] <- ifelse(is.na(x), NA_real_, vals)
    }
    # Drop the last column for restricted splines
    spline_terms <- spline_terms[, seq_len(n_knots - 1L), drop = FALSE]
  }

  # Normalize and return
  spline_terms / divisor
}


#' Build an additive design matrix for GAMs
#'
#' Constructs the expanded design matrix for generalized additive models by
#' appending spline basis terms according to per-column specifications. Each
#' column in `X` keeps its linear term; columns with a non-`NULL` specification
#' get additional spline basis columns appended.
#'
#' @param X Numeric matrix (n-by-b) of input covariates.
#' @param specifications A list of length `b` (number of columns in `X`).
#'   Each element is either `NULL` (no spline for that column) or a list with:
#'   \describe{
#'     \item{knots}{Numeric vector of knot locations (required).}
#'     \item{natural}{Logical, generate restricted splines? Default `TRUE`.}
#'     \item{power}{Numeric power for spline. Default `3`.}
#'     \item{penalty}{Numeric penalty for spline terms. Default `0`.}
#'     \item{normalized}{Logical, normalize spline terms? Default `FALSE`.}
#'   }
#' @param return_penalty Logical. If `TRUE`, return a list with both the
#'   design matrix and the penalty vector. Default `FALSE`.
#'
#' @returns If `return_penalty = FALSE`, a numeric matrix. If `TRUE`, a list
#'   with elements `X` (the design matrix) and `penalty` (numeric vector).
#'
#' @export
additive_design_matrix <- function(X, specifications, return_penalty = FALSE) {
  X <- as.matrix(X)
  n_obs <- nrow(X)
  n_cols <- ncol(X)

  # Handle single specification or NULL

  if (is.null(specifications)) {
    specifications <- vector("list", n_cols)
  } else if (is.list(specifications) && !is.null(specifications[["knots"]])) {
    # A single specification dict was provided; replicate for all columns
    specifications <- rep(list(specifications), n_cols)
  }

  if (length(specifications) != n_cols) {
    cli::cli_abort(
      "The number of specifications ({length(specifications)}) must match
       the number of columns in {.arg X} ({n_cols})."
    )
  }

  # Build columns and penalty vector

  result_cols <- list()
  penalties <- numeric(0)

  for (col_id in seq_len(n_cols)) {
    xvar <- X[, col_id]
    xspec <- specifications[[col_id]]

    # Always include the linear term
    result_cols <- c(result_cols, list(matrix(xvar, ncol = 1)))
    penalties <- c(penalties, 0)

    # Generate spline terms if specification is provided
    if (!is.null(xspec)) {
      spec <- process_spline_spec(xspec)
      spline_mat <- deli_spline(
        x = xvar,
        knots = spec$knots,
        power = spec$power,
        restricted = spec$natural,
        normalized = spec$normalized
      )
      result_cols <- c(result_cols, list(spline_mat))
      penalties <- c(penalties, rep(spec$penalty, ncol(spline_mat)))
    }
  }

  # Combine into a single matrix
  Xa <- do.call(cbind, result_cols)

  if (return_penalty) {
    list(X = Xa, penalty = penalties)
  } else {
    Xa
  }
}


#' Process and validate spline specification
#' @noRd
process_spline_spec <- function(spec) {
  if (is.null(spec[["knots"]])) {
    cli::cli_abort(
      "{.val knots} must be specified in the spline specification."
    )
  }

  # Fill in defaults for missing keys
  defaults <- list(
    natural = TRUE,
    power = 3,
    penalty = 0,
    normalized = FALSE
  )

  for (key in names(defaults)) {
    if (is.null(spec[[key]])) {
      spec[[key]] <- defaults[[key]]
    }
  }

  # Warn about unexpected keys
  expected_keys <- c("knots", "natural", "power", "penalty", "normalized")
  extra_keys <- setdiff(names(spec), expected_keys)
  if (length(extra_keys) > 0L) {
    cli::cli_warn(
      "Unexpected keys in spline specification: {.val {extra_keys}}.
       These will be ignored."
    )
  }

  spec
}

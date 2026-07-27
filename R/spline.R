#' Generate polynomial spline basis terms
#'
#' @description
#' Generates polynomial spline terms for a numeric vector at pre-specified
#' knot locations. Default is restricted (natural) cubic splines, but
#' unrestricted splines with different polynomial terms can also be generated.
#'
#' This function mirrors `spline()` in Python delicatessen, so code translated
#' from Python can keep its shape. Despite the name it has nothing to do with
#' [stats::spline()], which interpolates a curve through data points and returns
#' the interpolated values. `deli_spline()` builds a truncated power basis, a
#' matrix of design columns to be used as covariates in a regression. The `deli_`
#' prefix is there so that the two names do not collide.
#'
#' No base R function returns these truncated power columns, but the `splines`
#' package, which installs with every copy of R, builds the same spline spaces
#' in a B-spline parameterization. [splines::bs()] is the counterpart to
#' `restricted = FALSE`: for knots `k`,
#' `cbind(1, x, x^2, x^3, deli_spline(x, k, restricted = FALSE))` and
#' `cbind(1, splines::bs(x, knots = k))` span the same space and give identical
#' fitted values. [splines::ns()] is the nearest counterpart to
#' `restricted = TRUE` but not an exact one, since it also constrains the fit to
#' be linear beyond the boundary knots and so spans a subspace of the columns
#' here. `deli_spline()` is what deli offers for parity with Python
#' delicatessen, for the truncated power form itself, and because
#' [additive_design_matrix()] builds on it.
#'
#' @details
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
#' @param normalized Logical. If `TRUE`, divide spline terms by the range of
#'   knots (largest minus smallest) raised to `power`. With a single knot the
#'   range is zero, so the divisor is that knot raised to `power` instead.
#'   Default `FALSE`.
#'
#' @returns A matrix with `length(x)` rows. Number of columns is
#'   `length(knots)` for unrestricted or `length(knots) - 1` for restricted.
#'
#' @examples
#' # Restricted quadratic splines at four knots, so three basis columns
#' s <- deli_spline(1:59, knots = c(10, 20, 30, 40), power = 2)
#' dim(s)
#'
#' # Each term stays at zero below its knot and grows above it
#' s[c(5, 15, 25, 35, 45, 55), ]
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
#' @description
#' Constructs the expanded design matrix for generalized additive models by
#' appending spline basis terms according to per-column specifications. Each
#' column in `X` keeps its linear term; columns with a non-`NULL` specification
#' get additional spline basis columns appended.
#'
#' This function mirrors `additive_design_matrix()` in Python delicatessen, so
#' code translated from Python can keep its shape. There is no base R equivalent
#' of this construction, so this is the interface for it in deli as well.
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
#' @examples
#' set.seed(42)
#' X <- cbind(rnorm(50), rnorm(50))
#'
#' # The first column stays linear; the second also gets penalized splines
#' specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))
#' out <- additive_design_matrix(X, specs, return_penalty = TRUE)
#'
#' # Linear terms are unpenalized, spline terms carry the requested penalty
#' out$penalty
#'
#' dim(out$X)
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

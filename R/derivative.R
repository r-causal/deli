#' Numerical differentiation via finite differences
#'
#' Computes the Jacobian matrix of a vector-valued function using forward,
#' backward, or central difference approximation.
#'
#' @param func A function that takes a numeric vector and returns a numeric
#'   vector.
#' @param theta Numeric vector of parameter values at which to evaluate the
#'   derivative.
#' @param method Character string specifying the approximation method. One of
#'   `"capprox"` (central difference, default), `"fapprox"` (forward
#'   difference), or `"bapprox"` (backward difference).
#' @param dx Numeric step size for the finite difference (default `1e-9`).
#'
#' @return A matrix where element `[i, j]` is the partial derivative of the
#'   `i`-th output with respect to the `j`-th parameter.
#'
#' @keywords internal
approx_differentiation <- function(func, theta, method = "capprox", dx = 1e-9) {
  theta <- as.numeric(theta)
  p <- length(theta)
  # Identity matrix scaled by dx for perturbations
  shift <- diag(p) * dx

  # Evaluate func at each shifted parameter vector, return as matrix
  # Each row j is func(theta_shifted_j)
  generate_matrix <- function(x_shift) {
    # Drop dimensions on each evaluation so a column matrix (which any func
    # ending in %*% returns) becomes a length-n vector before stacking.
    rows <- lapply(seq_len(p), function(j) as.vector(func(x_shift[j, ])))
    do.call(rbind, rows)
  }

  if (method == "capprox") {
    # Central difference: (f(x+dx) - f(x-dx)) / (2*dx)
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f0 <- generate_matrix(lower)
    f1 <- generate_matrix(upper)
    deriv <- t(f1 - f0) / (2 * dx)

  } else if (method == "fapprox") {
    # Forward difference: (f(x+dx) - f(x)) / dx
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f_eval <- as.vector(func(theta))
    f0 <- matrix(rep(f_eval, each = p), nrow = p)
    f1 <- generate_matrix(upper)
    deriv <- t(f1 - f0) / dx

  } else if (method == "bapprox") {
    # Backward difference: (f(x) - f(x-dx)) / dx
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    f_eval <- as.vector(func(theta))
    f1 <- matrix(rep(f_eval, each = p), nrow = p)
    f0 <- generate_matrix(lower)
    deriv <- t(f1 - f0) / dx

  } else {
    cli::cli_abort("The method {.val {method}} is not supported.")
  }

  deriv
}

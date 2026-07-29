#' Numerical differentiation via finite differences
#'
#' Computes the Jacobian matrix of a vector-valued function using forward,
#' backward, or central difference approximation.
#'
#' @details
#' `dx` is an absolute perturbation, matching the `epsilon` argument of the
#' Python `delicatessen` library, so a `dx` carried over from there means the
#' same thing here. An absolute step has a limit, though: the spacing of the
#' doubles surrounding `theta` grows with `theta`, so a fixed `dx` spans fewer
#' and fewer of them as a parameter grows. At the default `dx` it spans about
#' 17,600 of them at `|theta| = 450`, 69 at `|theta| = 1.3e5`, and barely one
#' by `|theta| = dx / .Machine$double.eps`, about `4.5e6`. The perturbation is
#' rounded away with them, at first losing significant digits and past roughly
#' `1.7e7` leaving `theta + dx` equal to `theta`, at which point every
#' difference is zero and the Jacobian collapses.
#'
#' Each parameter's step is therefore floored at that magnitude's floating-point
#' resolution: it always spans at least ten thousand representable values, so
#' the perturbation actually applied reproduces the one intended to about four
#' significant digits. The floor is `1e4 * .Machine$double.eps * |theta|`, and
#' `.Machine$double.eps * |theta|` runs from one spacing to two as `|theta|`
#' climbs from one power of two to the next, so the floor engages wherever `dx`
#' would span fewer than ten to twenty thousand of them: every `|theta|` above
#' `dx / (1e4 * .Machine$double.eps)`, about `450` at the default `dx`. Below
#' that magnitude the step is `dx` exactly and nothing about the result changes.
#' Above it the step applied is the floor rather than `dx`, and the derivative
#' moves with it: at `|theta| = 500`, where `dx` still spans 17,592 values and
#' is resolved to five significant digits, the step is `1.11e-9`.
#'
#' Where the floor does engage, the derivative is accurate to a few parts in ten
#' thousand rather than to the roughly `1e-7` a representable step reaches.
#' `deriv_method = "exact"` takes no step at all and is unaffected by parameter
#' magnitude, so it is the better choice for a badly scaled problem.
#'
#' @param func A function that takes a numeric vector and returns a numeric
#'   vector.
#' @param theta Numeric vector of parameter values at which to evaluate the
#'   derivative.
#' @param method Character string specifying the approximation method. One of
#'   `"capprox"` (central difference, default), `"fapprox"` (forward
#'   difference), or `"bapprox"` (backward difference).
#' @param dx Numeric step size for the finite difference (default `1e-9`). The
#'   step is absolute, floored at the floating-point resolution of each
#'   parameter.
#'
#' @return A matrix where element `[i, j]` is the partial derivative of the
#'   `i`-th output with respect to the `j`-th parameter.
#'
#' @keywords internal
approx_differentiation <- function(func, theta, method = "capprox", dx = 1e-9) {
  theta <- as.numeric(theta)
  p <- length(theta)
  # A step of fewer than this many units in the last place of a parameter no
  # longer reproduces the perturbation asked for: rounding `theta + dx` back to
  # a representable double costs at most half an interval, so spanning n of them
  # holds the realized step to within 1/(2n) of `dx`.
  min_step_ulps <- 1e4
  step_size <- pmax(dx, min_step_ulps * .Machine$double.eps * abs(theta))
  # Identity matrix scaled by each parameter's step, for the partials
  shift <- diag(p) * step_size

  # Evaluate func at each shifted parameter vector, return as matrix
  # Each row j is func(theta_shifted_j)
  generate_matrix <- function(x_shift) {
    # Drop dimensions on each evaluation so a column matrix (which any func
    # ending in %*% returns) becomes a length-n vector before stacking.
    rows <- lapply(seq_len(p), function(j) as.vector(func(x_shift[j, ])))
    do.call(rbind, rows)
  }

  # The differences are divided before the transpose, where row j still belongs
  # to parameter j, so that each row meets its own step.
  if (method == "capprox") {
    # Central difference: (f(x+h) - f(x-h)) / (2*h)
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f0 <- generate_matrix(lower)
    f1 <- generate_matrix(upper)
    deriv <- t((f1 - f0) / (2 * step_size))
  } else if (method == "fapprox") {
    # Forward difference: (f(x+h) - f(x)) / h
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f_eval <- as.vector(func(theta))
    f0 <- matrix(rep(f_eval, each = p), nrow = p)
    f1 <- generate_matrix(upper)
    deriv <- t((f1 - f0) / step_size)
  } else if (method == "bapprox") {
    # Backward difference: (f(x) - f(x-h)) / h
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    f_eval <- as.vector(func(theta))
    f1 <- matrix(rep(f_eval, each = p), nrow = p)
    f0 <- generate_matrix(lower)
    deriv <- t((f1 - f0) / step_size)
  } else {
    cli::cli_abort(c(
      "The method {.val {method}} is not supported.",
      "i" = "Supported finite-difference options: {.val capprox},
             {.val fapprox}, {.val bapprox}. For exact automatic
             differentiation set {.code deriv_method = \"exact\"}."
    ))
  }

  deriv
}

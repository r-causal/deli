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
#' The floor rescues a step lost against a large `theta`. A difference is lost a
#' second way that no reading of the step can see, because nothing about the step
#' is wrong: where the step is applied exactly but the values of `func` are
#' large, the change it produces falls below the spacing of the doubles holding
#' them, both evaluations round to the same double, and the quotient is exactly
#' zero. That is indistinguishable, in the returned value, from a function that
#' is genuinely flat. What separates the two is the significance of the
#' difference against the magnitude of the values it was taken between, which is
#' what `deli_finite_difference_lost` reports; see **Value**.
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
#' @returns A matrix where element `[i, j]` is the partial derivative of the
#'   `i`-th output with respect to the `j`-th parameter. A Jacobian holding an
#'   entry whose difference was lost against the magnitude of the values of
#'   `func` is returned all the same, with a warning carrying the class
#'   `deli_finite_difference_lost`. One warning covers the call however many
#'   entries were lost, and its wording is the same at every call, so a caller
#'   wrapping several calls in `without_repeated_warnings()` reports it once for
#'   the operation.
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

  # Each branch leaves the two evaluations it differenced in `f0` and `f1` and
  # the quantity it divides by in `divisor`, so that the quotient and the reading
  # of it below are written once for the three methods. The divisor is a vector
  # of one entry per parameter, and the matrices carry parameter j in row j, so
  # recycling it down the rows meets each row with its own step.
  if (method == "capprox") {
    # Central difference: (f(x+h) - f(x-h)) / (2*h)
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f0 <- generate_matrix(lower)
    f1 <- generate_matrix(upper)
    divisor <- 2 * step_size
  } else if (method == "fapprox") {
    # Forward difference: (f(x+h) - f(x)) / h
    upper <- sweep(shift, 2, theta, FUN = function(s, th) th + s)
    f_eval <- as.vector(func(theta))
    f0 <- matrix(rep(f_eval, each = p), nrow = p)
    f1 <- generate_matrix(upper)
    divisor <- step_size
  } else if (method == "bapprox") {
    # Backward difference: (f(x) - f(x-h)) / h
    lower <- sweep(shift, 2, theta, FUN = function(s, th) th - s)
    f_eval <- as.vector(func(theta))
    f1 <- matrix(rep(f_eval, each = p), nrow = p)
    f0 <- generate_matrix(lower)
    divisor <- step_size
  } else {
    cli::cli_abort(c(
      "The method {.val {method}} is not supported.",
      "i" = "Supported finite-difference options: {.val capprox},
             {.val fapprox}, {.val bapprox}. For exact automatic
             differentiation set {.code deriv_method = \"exact\"}."
    ))
  }

  # The quotient is formed before the transpose, where row j still belongs to
  # parameter j.
  quotient <- (f1 - f0) / divisor
  warn_lost_difference(quotient, f0, f1, divisor)
  t(quotient)
}

#' Report a difference lost against the magnitude of the function values
#'
#' Reads each entry of the quotient against the smallest one the values it was
#' taken between can resolve. Neighboring doubles at a value `v` are
#' `.Machine$double.eps * |v|` apart, so a difference below that spacing is the
#' rounding of the two evaluations rather than any change in the function, and
#' dividing that spacing by the step gives the derivative the rounding alone is
#' worth: the noise floor of the entry. An entry no larger than its own noise
#' floor carries none of the digits of the values it came from.
#'
#' That reading on its own would condemn every honest zero, because a function
#' that does not change has a quotient of zero and a noise floor above it
#' whatever its values are. What separates the two is the size of the noise floor
#' against the scale of the Jacobian being built. A flat function evaluated at
#' values of order one has a floor of about `1e-7` at the default step, which is
#' negligible beside the derivatives an estimating equation produces, while a
#' collapse of the kind this reports leaves a floor of hundreds or thousands
#' beside a Jacobian of zeros. The two are separated by nine orders of magnitude
#' on the package's own fits, so the criterion is that the floor exceed the
#' largest derivative anywhere in the Jacobian, itself floored at one so that a
#' Jacobian which is genuinely all zeros is measured against something rather
#' than against nothing.
#'
#' Non-finite entries are excluded from both sides. A quotient that is `NA` or
#' infinite has a different problem, reported elsewhere, and the value of a
#' difference that was never taken says nothing about the scale of the rest.
#'
#' @param quotient The finite-difference quotients, one row per parameter.
#' @param f0,f1 The two evaluations each quotient was taken between, in the same
#'   shape.
#' @param divisor The quantity each row was divided by, one entry per parameter.
#'
#' @returns Invisible `NULL`, called for the warning.
#' @noRd
warn_lost_difference <- function(quotient, f0, f1, divisor) {
  noise <- .Machine$double.eps * pmax(abs(f0), abs(f1)) / divisor
  size <- abs(quotient)
  resolved <- is.finite(size)
  jacobian_scale <- max(1, size[resolved])
  lost <- resolved & is.finite(noise) & size <= noise & noise > jacobian_scale
  if (!any(lost)) {
    return(invisible(NULL))
  }
  # No value reaches the message, so every call raising it says the same thing
  # and an operation that differentiates more than once reports it once.
  cli::cli_warn(
    c(
      "!" = "A finite difference was lost to the magnitude of the function
             values, so at least one derivative came back as zero rather than
             as a number.",
      "i" = "The step changes the function by less than the spacing of the
             doubles holding its values, so both evaluations round to the same
             double and their difference carries none of their digits.",
      "i" = "Set {.code deriv_method = \"exact\"}, which takes no step, or
             rescale the problem. A larger {.arg dx} resolves the values too, at
             the cost of a coarser derivative."
    ),
    class = "deli_finite_difference_lost"
  )
  invisible(NULL)
}

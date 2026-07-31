#' M-Estimator
#'
#' S7 class for M-estimation via solving estimating equations with empirical
#' sandwich variance estimation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of parameters and n is the number of observations. Row names on
#'   that matrix name the parameters when `init` has none and every parameter is
#'   labeled; see [estimate()].
#' @param init Numeric vector of initial parameter values for the root-finding
#'   algorithm. Names on it label the parameters and take precedence over the
#'   row names of `stacked_equations`.
#' @inheritParams m_estimate
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction.
#' @param summed_equations A function that takes a numeric vector `theta` and
#'   returns the length-p vector of row sums of `stacked_equations` at `theta`,
#'   or `NULL` (default) to derive those sums from the full p-by-n return.
#'
#'   A fit reduces the estimating functions to those sums everywhere but the
#'   meat. The solver is given the summed equations to find a root of, and the
#'   bread is their Jacobian, so a fit builds the whole p-by-n matrix once per
#'   solver evaluation and once or twice more per parameter, all of it for
#'   arithmetic that is linear in it. An estimating function whose sums have a
#'   closed form, such as the \eqn{X^T r} of a regression score, can supply them
#'   here and give both steps only what they use. The meat is unaffected either
#'   way: it needs the per-observation contributions and takes the one full
#'   evaluation it always took, as does the validation at the starting values.
#'
#'   Under `deriv_method = "exact"` the reduction is called with a
#'   tangent-carrying `theta`, so it must be written in operations that carry
#'   derivatives: `t(X) %*% r` does, and [base::crossprod()] does not. See
#'   [auto_differentiation()] for which operations carry a tangent and where.
#'
#'   Anything that is neither `NULL` nor a function is refused here, with an
#'   error carrying the class `deli_summed_equations_error`. So is a return at
#'   the estimated values that is not numeric or does not hold one value per
#'   estimating equation, which [estimate()] reads.
#'
#'   [m_estimate()] and [gmm_estimate()] do not offer it. Both build
#'   `stacked_equations` themselves from a formula and the equation named in
#'   `.ee`, so the estimating functions a reduction would have to match are ones
#'   the caller never writes; this is for a system assembled by hand.
#' @param check_summed_equations Logical. When `TRUE` (default) and
#'   `summed_equations` was supplied, its value at the estimated values is
#'   compared against the row sums of the one full evaluation the meat is built
#'   from, and a disagreement raises an error carrying the class
#'   `deli_summed_equations_disagree`. The comparison costs one call to
#'   `summed_equations` and one reduction of a matrix the fit already holds.
#'
#'   What it is for is a reduction that sums some other system. The solver drives
#'   whatever it is given to zero, so such a reduction sends the fit to that
#'   other system's root and leaves the bread the Jacobian of one system and the
#'   meat the cross-product of another. Set it to `FALSE` to skip the
#'   comparison, which leaves a fit that reports estimates nobody solved for and
#'   a covariance with the shape of one and no claim to be one.
#'
#'   A reduction that is a multiple of the right one is not what the comparison
#'   sees. It vanishes where the right one does, so the fit lands at the same
#'   estimates and both quantities are at rounding where they are compared; the
#'   bread comes back as that multiple of the right bread. See
#'   [compute_sandwich()], whose argument of the same name reads the same
#'   comparison at a point the caller supplies.
#'
#' @returns An `MEstimator` S7 object. Call [estimate()] to solve the
#'   estimating equations and compute the sandwich variance.
#'
#' @export
#' @examples
#' # Estimating equations for the mean
#' y <- c(1, 2, 3, 4, 5)
#' psi <- function(theta) {
#'   matrix(y - theta[1], nrow = 1)
#' }
#' m <- MEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#' coef(m)
#' summary(m)
#'
#' # The solver and the bread only ever need the sums of the estimating
#' # equations, so an equation whose sums have a closed form can supply them
#' # and leave the whole matrix to the meat
#' summed <- function(theta) sum(y) - length(y) * theta[1]
#' MEstimator(stacked_equations = psi, init = 0, summed_equations = summed) |>
#'   estimate() |>
#'   summary()
MEstimator <- new_class(
  "MEstimator",
  parent = deli_estimator,
  constructor = function(
    stacked_equations,
    init,
    subset = NULL,
    finite_correction = NULL,
    summed_equations = NULL,
    check_summed_equations = TRUE
  ) {
    check_estimator_init(init)
    check_finite_correction(finite_correction)
    check_estimator_subset(subset, length(init))
    check_summed_function(summed_equations)
    check_summed_switch(check_summed_equations)

    init_names <- names(init)
    init <- as.numeric(init)
    names(init) <- init_names

    # Sort subset if provided
    if (!is.null(subset)) {
      subset <- sort(as.integer(subset))
    }

    new_object(
      S7_object(),
      stacked_equations = stacked_equations,
      summed_equations = summed_equations,
      check_summed_equations = check_summed_equations,
      init = init,
      subset = subset,
      finite_correction = finite_correction,
      n_params = length(init)
    )
  }
)

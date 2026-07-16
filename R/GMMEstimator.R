#' GMM Estimator
#'
#' S7 class for Generalized Method of Moments (GMM) estimation via
#' minimization of estimating equations with empirical sandwich variance
#' estimation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of estimating equations and n is the number of observations.
#'   The number of equations p must be greater than or equal to the number of
#'   parameters (length of `init`).
#' @param init Numeric vector of initial parameter values for the minimization
#'   algorithm.
#' @param subset Integer vector of parameter indices to solve for, or `NULL`
#'   (default) to solve for all parameters.
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction.
#' @param overid_maxiter Integer maximum iterations for the two-step iterative
#'   procedure for over-identified problems. Default `10L`.
#' @param overid_tolerance Numeric tolerance for convergence of the two-step
#'   iterative procedure for over-identified problems. Default `1e-9`.
#'
#' @returns A `GMMEstimator` S7 object. Call [estimate()] to minimize the
#'   estimating equations and compute the sandwich variance.
#'
#' @export
#' @examples
#' # Estimating equations for the mean (just-identified)
#' psi <- function(theta) {
#'   y <- c(1, 2, 3, 4, 5)
#'   matrix(y - theta[1], nrow = 1)
#' }
#' g <- GMMEstimator(stacked_equations = psi, init = c(0))
GMMEstimator <- new_class(
  "GMMEstimator",
  parent = deli_estimator,
  properties = list(
    overid_maxiter = new_property(
      class = class_integer,
      default = 10L
    ),
    overid_tolerance = new_property(
      class = class_numeric,
      default = 1e-9
    ),
    weight_matrix = new_property(
      class = NULL | class_double,
      default = NULL
    )
  ),
  constructor = function(stacked_equations, init, subset = NULL,
                         finite_correction = NULL,
                         overid_maxiter = 10L, overid_tolerance = 1e-9) {
    check_estimator_init(init)
    check_finite_correction(finite_correction)
    check_estimator_subset(subset, length(init))
    check_overid_scalar(overid_maxiter, "overid_maxiter")
    check_overid_scalar(overid_tolerance, "overid_tolerance")

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
      init = init,
      subset = subset,
      finite_correction = finite_correction,
      overid_maxiter = as.integer(overid_maxiter),
      overid_tolerance = overid_tolerance,
      n_params = length(init)
    )
  }
)

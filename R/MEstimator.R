#' M-Estimator
#'
#' S7 class for M-estimation via solving estimating equations with empirical
#' sandwich variance estimation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of parameters and n is the number of observations.
#' @param init Numeric vector of initial parameter values for the root-finding
#'   algorithm.
#' @param subset Integer vector of parameter indices to solve for, or `NULL`
#'   (default) to solve for all parameters.
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction.
#'
#' @returns An `MEstimator` S7 object. Call [estimate()] to solve the
#'   estimating equations and compute the sandwich variance.
#'
#' @export
#' @examples
#' # Estimating equations for the mean
#' psi <- function(theta) {
#'   y <- c(1, 2, 3, 4, 5)
#'   matrix(y - theta[1], nrow = 1)
#' }
#' m <- MEstimator(stacked_equations = psi, init = c(0))
MEstimator <- new_class(
  "MEstimator",
  parent = deli_estimator,
  constructor = function(
    stacked_equations,
    init,
    subset = NULL,
    finite_correction = NULL
  ) {
    check_estimator_init(init)
    check_finite_correction(finite_correction)
    check_estimator_subset(subset, length(init))

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
      n_params = length(init)
    )
  }
)

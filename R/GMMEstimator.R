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
#' y <- c(1, 2, 3, 4, 5)
#' psi <- function(theta) {
#'   matrix(y - theta[1], nrow = 1)
#' }
#' g <- GMMEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#' coef(g)
#'
#' # Instrumental variables with two instruments for a single treatment
#' # effect: 2 estimating equations, 1 parameter, so the system is
#' # over-identified and `MEstimator()` cannot solve it.
#' set.seed(42)
#' n <- 200
#' Z1 <- rbinom(n, 1, 0.5)
#' Z2 <- rnorm(n)
#' U <- rnorm(n)
#' A <- 0.5 * Z1 + 0.3 * Z2 + U + rnorm(n)
#' Y <- 2 * A - U + rnorm(n)
#'
#' psi_iv <- function(theta) {
#'   rbind(
#'     Z1 * (Y - theta[1] * A),
#'     Z2 * (Y - theta[1] * A)
#'   )
#' }
#'
#' # The two-step weight matrix update needs more than the default 10
#' # iterations to converge here.
#' g2 <- GMMEstimator(
#'   stacked_equations = psi_iv,
#'   init = 0,
#'   overid_maxiter = 200L
#' ) |>
#'   estimate()
#' summary(g2)
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
  constructor = function(
    stacked_equations,
    init,
    subset = NULL,
    finite_correction = NULL,
    overid_maxiter = 10L,
    overid_tolerance = 1e-9
  ) {
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

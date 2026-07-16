#' One-step M-estimation
#'
#' Creates an `MEstimator` and estimates it in one call, analogous to how
#' [stats::lm()] creates and fits a model in a single step. Supports both a
#' formula interface (for regression-family estimating equations) and a
#' function interface (for custom estimating equations).
#'
#' @param stacked_equations A formula or a function. When a formula, `data` and
#'   `.ee` must also be provided. When a function, it should take a numeric
#'   vector `theta` and return a p-by-n matrix.
#' @param data A data frame (required when `stacked_equations` is a formula).
#' @param .ee An estimating equation function that accepts `theta`, `X`, and
#'   the response as its third argument, plus optionally additional arguments
#'   (required when `stacked_equations` is a formula). The formula response is
#'   passed positionally, so it reaches whatever the function calls that
#'   argument (`y` for [ee_regression] or [ee_glm], `time` for [ee_aft]).
#' @param ... Additional arguments passed to `.ee`. When using the formula
#'   interface, these are evaluated with tidy evaluation in the context of
#'   `data`, so column names can be used directly (e.g., `event = status`).
#' @param init Numeric vector of initial parameter values. When `NULL`
#'   (default) and using the formula interface, a zero vector with names from
#'   the model matrix columns is generated automatically.
#' @param subset Integer vector of parameter indices to solve for, or `NULL`
#'   (default) to solve for all parameters. Indices are 1-based; parameters not
#'   listed are held fixed at their `init` values while the rest are solved.
#'   The variance estimator ignores `subset`. Passed straight through to the
#'   estimator constructor.
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction. When set, the meat
#'   matrix is rescaled and inference switches to the t-distribution with
#'   `df = n_obs - n_params`. Passed straight through to the estimator
#'   constructor.
#' @param solver Character string or function for the solver. Default `NULL`
#'   uses `"rootSolve"`.
#' @param maxiter Integer maximum iterations (default 5000).
#' @param tolerance Numeric convergence tolerance (default 1e-9).
#' @param deriv_method Character string for numerical differentiation method
#'   (default `"capprox"`).
#' @param dx Numeric step size for differentiation (default 1e-9).
#' @param allow_pinv Logical. Use pseudo-inverse if bread is singular? Default
#'   `TRUE`.
#'
#' @returns A fitted `MEstimator` object with populated `theta`, `variance`,
#'   etc. Use [coef()], [vcov()], [confint()], [summary()], or
#'   [broom::tidy()] to extract results.
#'
#' @export
#' @examples
#' # Formula interface
#' m <- m_estimate(mpg ~ wt + hp, data = mtcars,
#'                 .ee = ee_regression, model = "linear")
#' coef(m)
#' summary(m)
#'
#' # Function interface
#' y <- c(1, 2, 3, 4, 5)
#' m2 <- m_estimate(
#'   stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
#'   init = c(mean = 0)
#' )
#' coef(m2)
m_estimate <- function(stacked_equations, ...) {
  UseMethod("m_estimate")
}

#' @rdname m_estimate
#' @export
m_estimate.formula <- function(stacked_equations, data, .ee, ...,
                               init = NULL, subset = NULL,
                               finite_correction = NULL,
                               solver = NULL, maxiter = 5000,
                               tolerance = 1e-9, deriv_method = "capprox",
                               dx = 1e-9, allow_pinv = TRUE) {
  formula <- stacked_equations

  # Build model frame and extract X, y
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)

  # Evaluate ... with tidy evaluation in data context
  dots <- rlang::enquos(...)
  ee_args <- lapply(dots, function(q) rlang::eval_tidy(q, data = data))

  # Auto-generate init from model matrix if not provided
  if (is.null(init)) {
    init <- stats::setNames(rep(0, ncol(X)), colnames(X))
  }

  # Create the stacked_equations closure. The response is passed positionally
  # so it fills the first argument not named theta or X (`y`, `time`, ...).
  psi <- function(theta) {
    do.call(.ee, c(list(theta = theta, X = X), list(y), ee_args))
  }

  obj <- MEstimator(stacked_equations = psi, init = init, subset = subset,
                    finite_correction = finite_correction)
  estimate(obj, solver = solver, maxiter = maxiter,
           tolerance = tolerance, deriv_method = deriv_method,
           dx = dx, allow_pinv = allow_pinv)
}

#' @rdname m_estimate
#' @export
m_estimate.default <- function(stacked_equations, ..., init, subset = NULL,
                               finite_correction = NULL,
                               solver = NULL, maxiter = 5000,
                               tolerance = 1e-9, deriv_method = "capprox",
                               dx = 1e-9, allow_pinv = TRUE) {
  obj <- MEstimator(stacked_equations = stacked_equations, init = init,
                    subset = subset, finite_correction = finite_correction)
  estimate(obj, solver = solver, maxiter = maxiter,
           tolerance = tolerance, deriv_method = deriv_method,
           dx = dx, allow_pinv = allow_pinv)
}

#' One-step GMM estimation
#'
#' Creates a `GMMEstimator` and estimates it in one call. Supports both a
#' formula interface and a function interface, parallel to [m_estimate()].
#'
#' @inheritParams m_estimate
#' @param overid_maxiter Integer maximum iterations for the two-step iterative
#'   procedure for over-identified problems. Default `10L`.
#' @param overid_tolerance Numeric tolerance for convergence of the two-step
#'   iterative procedure. Default `1e-9`.
#'
#' @returns A fitted `GMMEstimator` object.
#'
#' @export
#' @examples
#' # Function interface
#' y <- c(1, 2, 3, 4, 5)
#' g <- gmm_estimate(
#'   stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
#'   init = c(mean = 0)
#' )
#' coef(g)
gmm_estimate <- function(stacked_equations, ...) {
  UseMethod("gmm_estimate")
}

#' @rdname gmm_estimate
#' @export
gmm_estimate.formula <- function(stacked_equations, data, .ee, ...,
                                 init = NULL, subset = NULL,
                                 finite_correction = NULL,
                                 solver = NULL, maxiter = 5000,
                                 tolerance = 1e-9, deriv_method = "capprox",
                                 dx = 1e-9, allow_pinv = TRUE,
                                 overid_maxiter = 10L,
                                 overid_tolerance = 1e-9) {
  formula <- stacked_equations

  # Build model frame and extract X, y
  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)

  # Evaluate ... with tidy evaluation in data context
  dots <- rlang::enquos(...)
  ee_args <- lapply(dots, function(q) rlang::eval_tidy(q, data = data))

  # Auto-generate init from model matrix if not provided
  if (is.null(init)) {
    init <- stats::setNames(rep(0, ncol(X)), colnames(X))
  }

  # Create the stacked_equations closure. The response is passed positionally
  # so it fills the first argument not named theta or X (`y`, `time`, ...).
  psi <- function(theta) {
    do.call(.ee, c(list(theta = theta, X = X), list(y), ee_args))
  }

  obj <- GMMEstimator(
    stacked_equations = psi, init = init, subset = subset,
    finite_correction = finite_correction,
    overid_maxiter = overid_maxiter, overid_tolerance = overid_tolerance
  )
  estimate(obj, solver = solver, maxiter = maxiter,
           tolerance = tolerance, deriv_method = deriv_method,
           dx = dx, allow_pinv = allow_pinv)
}

#' @rdname gmm_estimate
#' @export
gmm_estimate.default <- function(stacked_equations, ..., init, subset = NULL,
                                 finite_correction = NULL,
                                 solver = NULL, maxiter = 5000,
                                 tolerance = 1e-9, deriv_method = "capprox",
                                 dx = 1e-9, allow_pinv = TRUE,
                                 overid_maxiter = 10L,
                                 overid_tolerance = 1e-9) {
  obj <- GMMEstimator(
    stacked_equations = stacked_equations, init = init, subset = subset,
    finite_correction = finite_correction,
    overid_maxiter = overid_maxiter, overid_tolerance = overid_tolerance
  )
  estimate(obj, solver = solver, maxiter = maxiter,
           tolerance = tolerance, deriv_method = deriv_method,
           dx = dx, allow_pinv = allow_pinv)
}

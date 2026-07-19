#' Influence functions for M-Estimator
#'
#' Computes the influence function values for individual observations:
#' \deqn{IF(O_i; \theta) = B_n(\theta)^{-1} \psi(O_i; \theta)}
#'
#' @param object A fitted `MEstimator` object (after calling [estimate()]).
#' @param allow_pinv Logical. Use pseudo-inverse if bread is singular? Default
#'   `TRUE`.
#' @param ... Additional arguments (currently unused).
#'
#' @returns An n-by-p matrix of influence function values, where n is the
#'   number of observations and p is the number of parameters.
#'
#' @export
influence_functions <- new_generic(
  "influence_functions",
  "object",
  function(object, allow_pinv = TRUE, ...) {
    S7::S7_dispatch()
  }
)

method(influence_functions, deli_estimator) <- function(
  object,
  allow_pinv = TRUE,
  ...
) {
  check_estimated(object)

  # Evaluate estimating equations at theta-hat
  efunc_i <- object@stacked_equations(object@theta)
  if (is.null(dim(efunc_i))) {
    efunc_i <- matrix(efunc_i, nrow = 1)
  }

  # Invert bread matrix
  if (allow_pinv) {
    bread_inv <- tryCatch(
      solve(object@bread),
      error = function(e) {
        rlang::check_installed(
          "MASS",
          reason = "for pseudo-inverse when the bread matrix is singular."
        )
        MASS::ginv(object@bread)
      }
    )
  } else {
    bread_inv <- solve(object@bread)
  }

  # IF = B^{-1} * psi, return n-by-p
  t(bread_inv %*% efunc_i)
}

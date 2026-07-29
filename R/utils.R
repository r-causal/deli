#' Replace a `NULL` value with a default
#'
#' Returns `x` unless it is `NULL`, in which case it returns `y`.
#'
#' Defined here rather than taken from base R. `base::%||%` was added in R
#' 4.4.0, while DESCRIPTION declares a dependency on R (>= 4.3), so a package
#' that reaches for the base operator does not find it on every version it
#' claims to support. A binding in this namespace is found first on all of them
#' and removes the version dependency altogether.
#'
#' @param x A value to return when it is not `NULL`.
#' @param y The value to fall back on.
#' @returns `x`, or `y` when `x` is `NULL`.
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

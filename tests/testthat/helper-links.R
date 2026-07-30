# Helpers for pinning the derivative inverse_link() reports against its own
# inverse link. Two test files check the same contract from either side, one
# through predict() and one through ee_glm(), so the link roster and the central
# difference live here rather than being written out twice.

#' The canonical name of every link inverse_link() accepts
#'
#' The aliases (`"logistic"`, `"cauchy"`, `"square_root"`) are left out: each
#' resolves to a branch this roster already names.
#'
#' @return A character vector of link names.
inverse_link_names <- function() {
  c(
    "identity",
    "log",
    "logit",
    "probit",
    "cauchit",
    "loglog",
    "cloglog",
    "inverse",
    "sqrt"
  )
}

#' Differentiate an inverse link numerically by a central difference
#'
#' @param eta Numeric vector of linear predictor values.
#' @param link The link name, as accepted by `inverse_link()`.
#' @param step The half-width of the difference.
#'
#' @return A numeric vector of `d(mu)/d(eta)` values, one per element of `eta`.
numeric_inverse_link_deriv <- function(eta, link, step = 1e-6) {
  (inverse_link(eta + step, link)$mu - inverse_link(eta - step, link)$mu) /
    (2 * step)
}

#' Robust loss function derivatives
#'
#' @description
#' Computes the first derivative (psi function) of robust loss functions,
#' evaluated at the given residuals. Used internally by `ee_mean_robust()` and
#' `ee_robust_regression()`.
#'
#' This function mirrors `robust_loss_functions()` in Python delicatessen, so
#' code translated from Python can keep its shape. There is no base R equivalent
#' for these score functions, so this is the interface for them in deli as well.
#'
#' @param residuals Numeric vector of residuals.
#' @param loss Character string specifying the loss function. One of:
#'   `"huber"`, `"tukey"`, `"andrew"`, `"hampel"`, `"fair"`, `"cauchy"`,
#'   `"ullah"`, `"welsch"`.
#' @param k Numeric tuning constant. For `"hampel"`, a length-3 vector
#'   `c(a, b, c)` where `a < b < c`.
#'
#' @returns Numeric vector the same length as `residuals`.
#'
#' @export
#' @examples
#' r <- c(-5, -1, 0, 1, 5)
#' robust_loss_functions(r, "huber", k = 1.345)
#' robust_loss_functions(r, "tukey", k = 4.685)
robust_loss_functions <- function(residuals, loss, k) {
  if (!is.character(loss)) {
    cli::cli_abort("The {.arg loss} function should be a string.")
  }

  loss_l <- tolower(loss)

  if (loss_l == "huber") {
    # Clip residuals to [-k, k]
    pmin(pmax(residuals, -k), k)
  } else if (loss_l == "tukey") {
    ifelse(abs(residuals) <= k, residuals * (1 - (residuals / k)^2)^2, 0)
  } else if (loss_l == "andrew") {
    ifelse(abs(residuals) <= k * pi, sin(residuals / k), 0)
  } else if (loss_l == "fair") {
    residuals / (1 + abs(residuals) / k)
  } else if (loss_l == "cauchy") {
    residuals / (1 + (residuals / k)^2)
  } else if (loss_l == "ullah") {
    residuals * (1 + (residuals / k)^4)^(-2)
  } else if (loss_l == "welsch") {
    residuals * exp(-residuals^2 / (2 * k^2))
  } else if (loss_l == "hampel") {
    if (length(k) != 3) {
      cli::cli_abort(
        "The {.val hampel} loss function requires {.arg k} to be a
         length-3 vector {.code c(a, b, c)}."
      )
    }
    a <- k[1]
    b <- k[2]
    cc <- k[3]
    if (!(a < b && b < cc)) {
      cli::cli_abort(
        "The {.val hampel} loss function requires {.code a < b < c}."
      )
    }

    # Piecewise redescending score, written as nested conditionals so the
    # branch selection is tangent-aware under exact-mode autodiff (the masked
    # `ifelse` carries the selected arm's tangent). Each arm is evaluated over
    # the whole vector and selected by magnitude, which reproduces the earlier
    # logical-index assignment element for element.
    ar <- abs(residuals)
    sgn <- sign(residuals)
    ifelse(
      ar < a,
      # |r| < a: psi = r
      residuals,
      ifelse(
        ar < b,
        # a <= |r| < b: psi = a * sign(r)
        a * sgn,
        ifelse(
          ar < cc,
          # b <= |r| < c: psi = a * sign(r) * (c - |r|) / (c - b)
          a * sgn * (cc - ar) / (cc - b),
          # |r| >= c: psi = 0
          0
        )
      )
    )
  } else {
    cli::cli_abort("The loss function {.val {loss}} is not available.")
  }
}

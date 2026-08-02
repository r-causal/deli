# Estimating equation for regression calibration

Corrects for measurement error in a binary predictor using external
validation data. Scales the naive coefficient by the calibration factor.

## Usage

``` r
ee_regression_calibration(theta, beta, a, a_star, r, X = NULL, weights = NULL)
```

## Arguments

- theta:

  Numeric vector: the corrected coefficient first, then the calibration
  coefficients for the design `cbind(a_star, X)` (the `a_star`
  coefficient first). Length `2 + ncol(X)` when `X` is supplied, or
  length `3` when `X = NULL` (an intercept-only calibration model).

- beta:

  Numeric scalar. External estimate of the coefficient for the
  mismeasured predictor on the outcome.

- a:

  Numeric vector of gold-standard action measurements (validation sample
  only, where `r = 0`).

- a_star:

  Numeric vector of mismeasured action values.

- r:

  Numeric indicator: 0 for validation, 1 for main study.

- X:

  Optional design matrix for calibration model. Default `NULL` uses
  intercept only.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A `length(theta)`-by-n matrix (`2 + ncol(X)` rows, or `3` rows when
`X = NULL`). The first row is named `corrected_beta` and the second
`a_star`. The remaining rows are named `X_1` through `X_p` for the
columns of `X`, or `intercept` when `X = NULL`.

## Examples

``` r
# A binary exposure is measured with error in the main study. The external
# validation study regresses the gold-standard exposure on the mismeasured
# one, and that calibration slope rescales the naive outcome coefficient.
set.seed(789)
n <- 500
a_true <- rbinom(n, 1, 0.5)
a_star <- ifelse(a_true == 1, rbinom(n, 1, 0.85), rbinom(n, 1, 0.1))
r <- c(rep(0, 200), rep(1, 300))

# The gold standard is observed only in the validation sample, so the main
# study positions carry a 0 placeholder. It never reaches an estimate: the
# calibration model is multiplied by (1 - r).
a <- ifelse(r == 0, a_true, 0)

# `beta` is the naive coefficient for the mismeasured exposure, supplied here
# as a fixed external value. Stack an outcome model and pass its coefficient
# instead to propagate the uncertainty in that estimate as well.
psi <- function(theta) {
  ee_regression_calibration(theta, beta = 0.8, a = a, a_star = a_star, r = r)
}

m <- m_estimate(stacked_equations = psi, init = c(1, 0.1, 0.5))

# Corrected coefficient, then the calibration slope and intercept
coef(m)
#> corrected_beta         a_star      intercept 
#>      1.0853202      0.7371097      0.1764706 
```

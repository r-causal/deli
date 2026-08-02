# Estimating equations for weighted sensitivity analysis of the mean

Estimates the mean of an outcome with missing data using a weighted
sensitivity analysis approach. Handles MCAR, MAR, and MNAR mechanisms by
specifying a user-defined sensitivity function `q_eval` and a monotone
increasing distribution function `H_function`.

## Usage

``` r
ee_mean_sensitivity_analysis(theta, y, delta, X, q_eval, H_function)
```

## Arguments

- theta:

  Numeric vector of length `1 + b`, where `b` is the number of columns
  in `X`. The first element is the corrected mean; the remainder are
  regression coefficients.

- y:

  Numeric vector of n outcome values. Missing values should be indicated
  via the `delta` parameter.

- delta:

  Numeric vector of n indicators: 1 if `y` is observed, 0 if missing.
  Must not contain `NA`.

- X:

  Numeric n-by-b design matrix for the missingness model. Should include
  an intercept column. Must not contain `NA`.

- q_eval:

  Numeric vector of n evaluated sensitivity function values, i.e.
  `q(Y, alpha)`.

- H_function:

  A function mapping real values to `[0, 1]` that is monotone increasing
  (e.g.,
  [`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md)).

## Value

A `(1+b)`-by-n matrix of estimating equation contributions, with the
first row named `corrected_mean` and the missingness model rows named
`X_1` through `X_b` for the columns of `X`.

## Examples

``` r
# An outcome observed for only part of the sample, with missingness driven by
# the measured covariate W.
set.seed(42)
n <- 500
W <- rbinom(n, 1, 0.5)
Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
delta <- rbinom(n, 1, inverse_logit(2 + W))

# Missing outcomes never enter the estimating equation, so any placeholder
# value works; a zero keeps the arithmetic finite.
Y <- ifelse(delta == 1, Y_full, 0)
X <- cbind(1, W) # Missingness model design matrix

# A sensitivity function of zero everywhere assumes the outcome is missing at
# random given W. Nonzero values encode departures from that assumption.
psi <- function(theta) {
  ee_mean_sensitivity_analysis(
    theta,
    y = Y,
    delta = delta,
    X = X,
    q_eval = rep(0, n),
    H_function = inverse_logit
  )
}

# theta holds the corrected mean, started near the observed outcomes,
# followed by the two missingness model coefficients.
m <- m_estimate(stacked_equations = psi, init = c(180, 0, 0))
coef(m)
#> corrected_mean            X_1            X_2 
#>     183.174625       1.698074       1.061031 
```

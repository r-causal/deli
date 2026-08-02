# Generate predicted values from a regression model

Computes predicted outcomes, their variance, and Wald-type confidence
intervals from estimated regression coefficients and their covariance
matrix. This is a post-processing utility meant to be used after
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
has been fitted.

## Usage

``` r
regression_predictions(X, theta, covariance, offset = NULL, alpha = 0.05)
```

## Arguments

- X:

  Numeric n-by-p design matrix of covariate values for prediction.

- theta:

  Numeric vector of p estimated coefficients (from `coef(m)`).

- covariance:

  Numeric p-by-p covariance matrix (from `vcov(m)`).

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

- alpha:

  Numeric significance level for confidence intervals. Default `0.05`
  (95% CIs).

## Value

A data frame with n rows and columns: `predicted`, `variance`, `lower`,
`upper`. The rows are labeled by position regardless of the row names of
`X`, which say nothing about the predictions made from it.

## Details

No transformations are applied. For logistic models this returns
log-odds (not probabilities). Apply
[`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html) for the
probability scale, or
[`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md)
if the values feed a transform passed to
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md)
with `deriv_method = "exact"`.

## Examples

``` r
set.seed(1)
n <- 200
dat <- data.frame(x = rnorm(n), z = rbinom(n, 1, 0.5))
dat$y <- 1 + 0.5 * dat$x + 2 * dat$z + rnorm(n)

m <- m_estimate(y ~ x + z, data = dat, .ee = ee_regression, model = "linear")

# Predict along a small grid of x, holding z at 1. The columns of the grid
# must match the order of the coefficients, intercept first.
X_new <- cbind(1, x = c(-1, 0, 1), z = 1)
regression_predictions(X_new, theta = coef(m), covariance = vcov(m))
#>   predicted   variance    lower    upper
#> 1  2.501293 0.02222253 2.209117 2.793469
#> 2  3.065570 0.01192616 2.851528 3.279612
#> 3  3.629847 0.01613069 3.380919 3.878776
```

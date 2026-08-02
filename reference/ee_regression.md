# Estimating equation for regression

Returns a p-by-n matrix of estimating equation contributions for
regression models. Supports linear, logistic, and Poisson regression:
\$\$\psi_i(\theta) = \\Y_i - g(X_i^T \theta)\\ X_i\$\$

## Usage

``` r
ee_regression(theta, X, y, model, weights = NULL, offset = NULL)
```

## Arguments

- theta:

  Numeric vector of length p (number of covariates).

- X:

  Numeric n-by-p design matrix.

- y:

  Numeric vector of n observed outcome values.

- model:

  Character string: `"linear"`, `"logistic"`, or `"poisson"`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix.

## Examples

``` r
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_regression,
  model = "linear"
)
summary(fit)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 32
#> Parameters: 3
#> 
#>               Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> (Intercept)    37.2273     1.9389    19.2000    33.4271    41.0275     <2e-16   270.5102
#> wt             -3.8778     0.6199    -6.2553    -5.0929    -2.6628   3.97e-10    31.2310
#> hp             -0.0318     0.0066    -4.7807    -0.0448    -0.0187   1.75e-06    19.1270

# The same equation fits a logistic regression through the model argument.
fit_logit <- m_estimate(
  vs ~ mpg,
  data = mtcars,
  .ee = ee_regression,
  model = "logistic"
)
coef(fit_logit)
#> (Intercept)         mpg 
#>  -8.8330726   0.4304135 
```

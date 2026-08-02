# Estimating equation for bridge penalized regression

Returns a p-by-n matrix for bridge penalized regression. Bridge is the
general case: ridge is `gamma = 2`, LASSO approximation is
`gamma = 1 + epsilon`.

## Usage

``` r
ee_bridge_regression(
  theta,
  X,
  y,
  model,
  penalty,
  gamma,
  weights = NULL,
  center = 0,
  offset = NULL
)
```

## Arguments

- theta:

  Numeric vector of length p.

- X:

  Numeric n-by-p design matrix.

- y:

  Numeric vector of n observed outcome values.

- model:

  Character string: `"linear"`, `"logistic"`, or `"poisson"`.

- penalty:

  Numeric scalar or vector of length p. Must be non-negative.

- gamma:

  Numeric bridge exponent. Must be at least 1. Values below 2 yield a
  penalty that is not everywhere differentiable, so the sandwich
  variance is not defined in all settings and a warning is issued.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- center:

  Numeric scalar or vector. Default `0`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix.

## Examples

``` r
# A penalty vector gives one value per column of the design matrix. A scalar
# penalty would shrink the intercept along with the slopes. The estimating
# equation carries the penalty's derivative, which varies like
# |theta|^(gamma - 1) near the penalty center. That derivative is itself
# differentiable only once gamma reaches 2, so gamma = 2.3 issues no
# warning while a value below 2 would.
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_bridge_regression,
  model = "linear",
  penalty = c(0, 5, 5),
  gamma = 2.3
)
coef(fit)
#> (Intercept)          wt          hp 
#> 35.17711995 -2.76260075 -0.04225662 
```

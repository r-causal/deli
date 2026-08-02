# Estimating equation for differentiable LASSO regression

Uses a smooth approximation to the L1 penalty based on the standard
normal CDF and PDF.

## Usage

``` r
ee_dlasso_regression(
  theta,
  X,
  y,
  model,
  penalty,
  s = 1e-06,
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

- s:

  Numeric smoothing parameter. Must be greater than zero. Default
  `1e-6`.

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
# equation carries the penalty's derivative. Here that derivative is smooth,
# so the estimating equation is differentiable and no warning is issued,
# unlike ee_lasso_regression() whose derivative has unbounded slope at the
# penalty center.
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_dlasso_regression,
  model = "linear",
  penalty = c(0, 5, 5)
)
coef(fit)
#> (Intercept)          wt          hp 
#> 36.68027321 -3.58300419 -0.03451029 
```

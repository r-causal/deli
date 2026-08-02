# Estimating equation for approximate LASSO regression

Uses the bridge penalty with `gamma = 1 + epsilon` to approximate LASSO
(L1 penalty). The true LASSO is not differentiable at zero, so an
approximation is used.

## Usage

``` r
ee_lasso_regression(
  theta,
  X,
  y,
  model,
  penalty,
  epsilon = 0.003,
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

- epsilon:

  Numeric approximation parameter. Default `0.003`.

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
# penalty would shrink the intercept along with the slopes.
#
# The approximate L1 penalty enters the estimating equation as its own
# derivative, and that derivative has unbounded slope at the penalty center.
# The estimating equation is therefore not differentiable there, so the
# bread matrix is undefined and the fit warns once that the sandwich
# variance should not be trusted here.
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_lasso_regression,
  model = "linear",
  penalty = c(0, 5, 5)
)
#> Warning: The estimating equation for the chosen penalized regression model is not always
#> differentiable. Therefore, the bread matrix is not always defined for finite
#> samples, and the sandwich should not be used to estimate the variance.
coef(fit)
#> (Intercept)          wt          hp 
#> 36.67817397 -3.58183470 -0.03452163 
```

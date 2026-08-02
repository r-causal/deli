# Estimating equation for elastic net regression

Combines L1 (approximate LASSO via bridge) and L2 (ridge) penalties at a
given ratio. When `ratio = 1`, this is LASSO; when `ratio = 0`, ridge.

## Usage

``` r
ee_elasticnet_regression(
  theta,
  X,
  y,
  model,
  penalty,
  ratio,
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

- ratio:

  Numeric between 0 and 1. Proportion of L1 vs L2 penalty.

- epsilon:

  Numeric LASSO approximation parameter. Default `0.003`.

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
# The L1 half of the penalty enters the estimating equation as its own
# derivative, and that derivative has unbounded slope at the penalty center.
# The estimating equation is therefore not differentiable there, so the
# bread matrix is undefined and the fit warns once that the sandwich
# variance should not be trusted here.
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_elasticnet_regression,
  model = "linear",
  penalty = c(0, 5, 5),
  ratio = 0.5
)
#> Warning: The estimating equation for the chosen penalized regression model is not always
#> differentiable. Therefore, the bread matrix is not always defined for finite
#> samples, and the sandwich should not be used to estimate the variance.
coef(fit)
#> (Intercept)          wt          hp 
#> 36.06463913 -3.24677056 -0.03768788 
```

# Estimating equation for ridge regression

Returns a p-by-n matrix of estimating equation contributions for ridge
(L2-penalized) regression: \$\$\psi_i(\theta) = \\Y_i - g(X_i^T
\theta)\\ X_i - \frac{\lambda}{n} \theta\$\$

## Usage

``` r
ee_ridge_regression(
  theta,
  X,
  y,
  model,
  penalty,
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

  Numeric scalar or vector of length p. Must be non-negative. Penalty
  terms scaled by n internally.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- center:

  Numeric scalar or vector. Center for the penalty. Default `0`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix.

## Examples

``` r
# A penalty vector gives one value per column of the design matrix. A scalar
# penalty would shrink the intercept along with the slopes.
fit <- m_estimate(
  mpg ~ wt + hp,
  data = mtcars,
  .ee = ee_ridge_regression,
  model = "linear",
  penalty = c(0, 5, 5)
)
coef(fit)
#> (Intercept)          wt          hp 
#> 35.59224774 -2.98849397 -0.04013219 
```

# Estimating equation for robust regression

Returns a p-by-n matrix for robust regression using the specified loss
function (linear regression only): \$\$\psi_i(\theta) = f_k(Y_i - X_i^T
\theta) X_i\$\$

## Usage

``` r
ee_robust_regression(
  theta,
  X,
  y,
  model,
  k,
  loss = "huber",
  weights = NULL,
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

  Character string. Currently only `"linear"` is supported.

- k:

  Numeric tuning parameter for the loss function.

- loss:

  Character string specifying the loss function. Default `"huber"`. See
  [`robust_loss_functions()`](https://r-causal.github.io/deli/reference/robust_loss_functions.md).

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix.

## Examples

``` r
# The Huber loss is convex, so its estimating function has a single root and
# seeding the search is about reaching that root rather than choosing among
# several. What makes the seed necessary is that the Huber psi is bounded:
# far from the solution every residual is past the tuning constant k, every
# contribution saturates at k, and the estimating function is constant with a
# Jacobian of exactly zero. Starting from zero lands in that flat region,
# where the solver has no slope to follow, so start from a least-squares fit.
start <- coef(lm(weight ~ height, data = robust_regress))

fit <- m_estimate(
  weight ~ height,
  data = robust_regress,
  .ee = ee_robust_regression,
  model = "linear",
  k = 1.345,
  loss = "huber",
  init = start
)
coef(fit)
#> (Intercept)      height 
#>  -36.791469    0.619315 
```

# Robust loss function derivatives

Computes the first derivative (psi function) of robust loss functions,
evaluated at the given residuals. Used internally by
[`ee_mean_robust()`](https://r-causal.github.io/deli/reference/ee_mean_robust.md)
and
[`ee_robust_regression()`](https://r-causal.github.io/deli/reference/ee_robust_regression.md).

This function mirrors `robust_loss_functions()` in Python delicatessen,
so code translated from Python can keep its shape. There is no base R
equivalent for these score functions, so this is the interface for them
in deli as well.

## Usage

``` r
robust_loss_functions(residuals, loss, k)
```

## Arguments

- residuals:

  Numeric vector of residuals.

- loss:

  Character string specifying the loss function. One of: `"huber"`,
  `"tukey"`, `"andrew"`, `"hampel"`, `"fair"`, `"cauchy"`, `"ullah"`,
  `"welsch"`.

- k:

  Numeric tuning constant. For `"hampel"`, a length-3 vector
  `c(a, b, c)` where `a < b < c`.

## Value

Numeric vector the same length as `residuals`.

## Examples

``` r
r <- c(-5, -1, 0, 1, 5)
robust_loss_functions(r, "huber", k = 1.345)
#> [1] -1.345 -1.000  0.000  1.000  1.345
robust_loss_functions(r, "tukey", k = 4.685)
#> [1]  0.0000000 -0.9109563  0.0000000  0.9109563  0.0000000
```

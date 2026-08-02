# Compute confidence bands from theta and covariance

Compute confidence bands from theta and covariance

## Usage

``` r
compute_confidence_bands(
  theta,
  covariance,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL
)
```

## Arguments

- theta:

  Numeric parameter vector.

- covariance:

  Numeric covariance matrix.

- alpha:

  Significance level. Default `0.05`.

- method:

  `"supt"` or `"bonferroni"`. Default `"supt"`.

- n_draws:

  Number of MVN draws for sup-t. Default `1e5`. See
  [`confidence_bands()`](https://r-causal.github.io/deli/reference/confidence_bands.md)
  for why this differs from Python's `1e6` estimator default and how to
  match it.

- seed:

  RNG seed. Default `NULL`.

## Value

A p-by-2 matrix with columns `"lower"` and `"upper"`. Rows take their
names from `theta`, when it has any.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

# Bands from the estimates and covariance alone, without the fitted object
compute_confidence_bands(coef(fit), covariance = vcov(fit),
                         method = "supt", seed = 1)
#>                   lower       upper
#> (Intercept) 32.71991838 41.73462185
#> wt          -5.31895808 -2.43670340
#> hp          -0.04722284 -0.01632305
```

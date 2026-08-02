# Broom tidiers for deli estimators

[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
methods for `MEstimator` and `GMMEstimator` objects. These allow deli
results to flow into tidyverse pipelines. Requires the generics package
(or broom) to be loaded.

## Arguments

- x:

  A fitted `MEstimator` or `GMMEstimator` object.

- conf.int:

  Logical. Include confidence intervals? Default `FALSE`.

- conf.level:

  Numeric confidence level for intervals. Default `0.95`.

- ...:

  Not used. `tidy()` requires them to be empty, so that a misspelled
  `conf.int` or `conf.level` is an error rather than a table silently
  returned without intervals or at the default level. `glance()` has no
  optional argument for a wrong name to displace and ignores them.

## Value

- `tidy()`: A data.frame with columns `term`, `estimate`, `std.error`,
  `statistic`, `p.value`, `s.value`. If `conf.int = TRUE`, also includes
  `conf.low` and `conf.high`. A `p.value` that underflows to exactly
  zero is reported as `0` alongside an infinite `s.value`; see
  [`s_values()`](https://r-causal.github.io/deli/reference/s_values.md).

- `glance()`: A single-row data.frame with model-level summaries:
  `nobs`, `npar`, `estimator`, `finite_correction`, and the Hansen
  J-statistic of an over-identified GMM fit in `j_statistic`, `j_df` and
  `j_p_value`. The three J columns are present on every fit and hold the
  typed missing value of their own type where there is no such
  statistic, which is every M-estimation fit and every just-identified
  or `subset` GMM fit; see
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  for what the statistic reads and where it is left unset.

## See also

[deli-augment](https://r-causal.github.io/deli/reference/deli-augment.md),
the third broom generic, which returns the observation-level fitted
values, intervals, and residuals.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

# Attaching generics or broom makes the shorter tidy(fit) and glance(fit)
# calls work as well
generics::tidy(fit, conf.int = TRUE)
#>          term    estimate   std.error statistic      p.value   s.value
#> 1 (Intercept) 37.22727012 1.938920903 19.199994 3.701332e-82 270.51006
#> 2          wt -3.87783074 0.619927640 -6.255296 3.967632e-10  31.23100
#> 3          hp -0.03177295 0.006646058 -4.780721 1.746680e-06  19.12695
#>      conf.low   conf.high
#> 1 33.42705498 41.02748525
#> 2 -5.09286659 -2.66279490
#> 3 -0.04479898 -0.01874691

generics::glance(fit)
#>   nobs npar  estimator finite_correction j_statistic j_df j_p_value
#> 1   32    3 MEstimator              <NA>          NA   NA        NA
```

# Z-scores for M-Estimator parameters

Computes Wald-type Z-scores: \\(\hat{\theta} - \theta_0) /
\widehat{SE}(\hat{\theta})\\.

This function mirrors `m.z_scores()` in Python delicatessen, so code
translated from Python can keep its shape.

## Usage

``` r
z_scores(object, null = 0, ...)
```

## Arguments

- object:

  A fitted `MEstimator` object (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- null:

  Numeric null hypothesis value(s). Default `0`.

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

A numeric vector of Z-scores.

## See also

[`summary()`](https://r-causal.github.io/deli/reference/deli-display.md)
and
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
which report the same Z-scores in table form alongside the other
parameter-level results.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

z_scores(fit)
#> (Intercept)          wt          hp 
#>   19.199998   -6.255299   -4.780721 
```

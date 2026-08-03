# S-values (surprisal) for M-Estimator parameters

Computes Shannon Information values (S-values) as \\S = -\log_2(P)\\,
where \\P\\ is the corresponding P-value.

This function mirrors `m.s_values()` in Python delicatessen, so code
translated from Python can keep its shape.

## Usage

``` r
s_values(object, null = 0, ...)
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

A numeric vector of S-values.

## Details

A P-value small enough to underflow to exactly zero has an S-value of
`Inf`. That is the limit the surprisal is heading toward rather than a
defect. The smallest P a double can hold is \\2^{-1074}\\, so a P that
arrives as zero stands for more than a thousand bits of surprisal, past
the range a double can name. The infinity reports evidence beyond
measurement, where any finite substitute would name a number the fit
does not support.

## See also

[`summary()`](https://r-causal.github.io/deli/reference/deli-display.md)
and
[`tidy()`](https://r-causal.github.io/deli/reference/deli-tidiers.md),
which report the same S-values in table form alongside the other
parameter-level results.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

s_values(fit)
#> (Intercept)          wt          hp 
#>   270.51017    31.23104    19.12696 
```

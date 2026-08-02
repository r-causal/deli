# P-values for M-Estimator parameters

Computes two-sided Wald-type P-values from Z-scores. The Z-scores are
compared to the standard normal distribution by default. When a
`finite_correction` is set on the fit, they are compared instead to the
t-distribution with \\n - p\\ degrees of freedom.

This function mirrors `m.p_values()` in Python delicatessen, so code
translated from Python can keep its shape.

## Usage

``` r
p_values(object, null = 0, ...)
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

A numeric vector of P-values.

## See also

[`summary()`](https://r-causal.github.io/deli/reference/deli-display.md)
and
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html),
which report the same P-values in table form alongside the other
parameter-level results.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

p_values(fit)
#>  (Intercept)           wt           hp 
#> 3.701056e-82 3.967539e-10 1.746674e-06 
```

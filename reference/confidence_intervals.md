# Confidence intervals for M-Estimator parameters

Computes two-sided Wald-type \\(1 - \alpha) \times 100\\\\ confidence
intervals using the point estimates and sandwich variance:
\\\hat{\theta} \pm c\_{\alpha/2} \times \widehat{SE}(\hat{\theta})\\.
The critical value \\c\_{\alpha/2}\\ comes from the standard normal
distribution by default. When a `finite_correction` is set on the fit,
it comes instead from the t-distribution with \\n - p\\ degrees of
freedom, matching the finite-sample adjustment of the variance.

This function mirrors `m.confidence_intervals()` in Python delicatessen,
so code translated from Python can keep its shape.

## Usage

``` r
confidence_intervals(object, alpha = 0.05, ...)
```

## Arguments

- object:

  A fitted `MEstimator` object (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- alpha:

  Numeric significance level, between 0 and 1. Default `0.05` for 95%
  confidence intervals.

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

A p-by-2 matrix with columns `"lower"` and `"upper"`.

## See also

[`confint()`](https://r-causal.github.io/deli/reference/deli-generics.md),
the standard R accessor for the same intervals. It returns identical
values but is parameterized by `level = 0.95` where this function takes
`alpha = 0.05`.

## Examples

``` r
psi <- function(theta) {
  y <- c(1, 2, 3, 4, 5)
  matrix(y - theta[1], nrow = 1)
}
m <- m_estimate(stacked_equations = psi, init = 0)
confidence_intervals(m)
#>           lower   upper
#> theta_1 1.76041 4.23959
```

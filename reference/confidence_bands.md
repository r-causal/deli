# Confidence bands for parameter vectors

Computes simultaneous confidence bands that provide coverage for the
entire parameter vector, adjusting for multiple comparisons. The formula
is: \$\$\hat{\theta} \pm \hat{c}\_{\alpha/2} \times
\widehat{SE}(\hat{\theta})\$\$ where \\\hat{c}\\ is the adjusted
critical value.

## Usage

``` r
confidence_bands(
  object,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL,
  subset = NULL,
  covariance = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `MEstimator` object, or a numeric vector of parameter
  estimates.

- alpha:

  Numeric significance level, between 0 and 1. Default `0.05`.

- method:

  Character string. `"supt"` (default) for supremum-t or `"bonferroni"`
  for Bonferroni correction.

- n_draws:

  Integer number of MVN draws for the sup-t method. Default `1e5`. The
  critical value is a Monte Carlo estimate of a fixed sup-t quantile,
  and at `1e5` draws the band half-width varies by roughly 0.2% across
  independent draws, which is negligible relative to the band width.
  Python's estimator method defaults to `1e6`; both defaults target the
  same quantity and agree to within Monte Carlo error, so `deli` keeps
  the smaller default for faster computation. Set `n_draws = 1e6` to
  match Python's estimator default. Note that seeded band values are not
  reproducible across the two languages because the MVN samplers differ.

- seed:

  Integer seed for reproducibility. Default `NULL`.

- subset:

  Integer vector of parameter indices to compute bands for. Default
  `NULL` (all parameters).

- covariance:

  Numeric covariance matrix (only when `object` is numeric).

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

A p-by-2 matrix with columns `"lower"` and `"upper"`. For a fitted
estimator the rows are named for the parameters, as in
[`confint()`](https://r-causal.github.io/deli/reference/deli-generics.md).
For a numeric `object` they take their names from it, when it has any.

## Examples

``` r
# Two independent samples, each contributing one mean parameter
set.seed(42)
n <- 200
y1 <- rnorm(n, 2)
y2 <- rnorm(n, 3)

psi <- function(theta) {
  rbind(y1 - theta[1], y2 - theta[2])
}

m <- m_estimate(stacked_equations = psi, init = c(0, 0))

# Simultaneous bands cover both means at once, so they are wider than the
# pointwise intervals from confint(m)
confidence_bands(m, method = "supt", seed = 1)
#>            lower    upper
#> theta_1 1.818388 2.126643
#> theta_2 2.861506 3.161062
```

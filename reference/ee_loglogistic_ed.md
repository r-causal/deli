# Estimating equation for delta-effective dose (log-logistic)

Computes the effective dose at level delta from the log-logistic model.
Should be stacked with
[`ee_loglogistic()`](https://r-causal.github.io/deli/reference/ee_loglogistic.md).

## Usage

``` r
ee_loglogistic_ed(theta, dose, delta, lower, upper, ed50, steepness)
```

## Arguments

- theta:

  Numeric scalar. The ED(delta) parameter.

- dose:

  Numeric vector (used for dimension only).

- delta:

  Numeric effective dose level of interest.

- lower:

  Numeric lower limit parameter.

- upper:

  Numeric upper limit parameter.

- ed50:

  Numeric estimated ED50.

- steepness:

  Numeric steepness parameter.

## Value

A 1-by-n matrix.

## See also

[`ee_loglogistic()`](https://r-causal.github.io/deli/reference/ee_loglogistic.md),
the dose-response equation this one is stacked with.

## Examples

``` r
# The lower limit is held at zero instead of being estimated: root length
# cannot fall below zero, and the five-parameter stack diverges on these
# data. The lower-limit row of the log-logistic equation is therefore
# dropped, leaving the upper limit, ED50, and steepness to be estimated
# alongside the ED90.
psi <- function(theta) {
  loglogistic <- ee_loglogistic(
    c(0, theta[1:3]),
    dose = inderjit$dose,
    response = inderjit$response
  )
  ed90 <- ee_loglogistic_ed(
    theta[4],
    dose = inderjit$dose,
    delta = 0.9,
    lower = 0,
    upper = theta[1],
    ed50 = theta[2],
    steepness = theta[3]
  )
  rbind(loglogistic[-1, , drop = FALSE], ed90)
}

# This stacked system is sensitive to its starting values, so a reasonable
# init matters more here than it does for most equations.
m <- m_estimate(stacked_equations = psi, init = c(8, 2, 1, 5))

# Upper limit, ED50, steepness, and the ED90. Stacking is what gives the ED90
# a standard error of its own, so summary() rather than coef() shows the
# gain.
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 24
#> Parameters: 4
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1     7.8554     0.1540    51.0194     7.5537     8.1572     <2e-16        Inf
#> theta_2     3.2634     0.2657    12.2813     2.7426     3.7842     <2e-16   112.7545
#> theta_3     2.4703     0.2924     8.4490     1.8973     3.0434     <2e-16    54.9177
#> theta_4     7.9423     1.2765     6.2220     5.4405    10.4442   4.91e-10    30.9244
```

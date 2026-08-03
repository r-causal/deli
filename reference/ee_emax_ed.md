# Estimating equation for delta-effective dose (E-max)

Computes the effective dose at level delta from the E-max model. Should
be stacked with
[`ee_emax()`](https://r-causal.github.io/deli/reference/ee_emax.md).

## Usage

``` r
ee_emax_ed(theta, dose, delta, ed50)
```

## Arguments

- theta:

  Numeric scalar. The ED(delta) parameter.

- dose:

  Numeric vector (used for dimension only).

- delta:

  Numeric effective dose level of interest.

- ed50:

  Numeric estimated ED50 from E-max model.

## Value

A 1-by-n matrix.

## See also

[`ee_emax()`](https://r-causal.github.io/deli/reference/ee_emax.md), the
dose-response equation this one is stacked with.

## Examples

``` r
# This equation carries no information on its own, so stack it with the
# E-max equation that supplies the ED50. Here delta = 0.9 requests the ED90.
psi <- function(theta) {
  emax <- ee_emax(
    theta[1:3],
    dose = inderjit$dose,
    response = inderjit$response
  )
  ed90 <- ee_emax_ed(
    theta[4],
    dose = inderjit$dose,
    delta = 0.9,
    ed50 = theta[3]
  )
  rbind(emax, ed90)
}

m <- m_estimate(stacked_equations = psi, init = c(8, -8, 2, 10))

# theta_4 is the ED90. Stacking is what gives it a standard error of its own,
# so summary() rather than coef() is what shows the gain.
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 24
#> Parameters: 4
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1     8.2151     0.2540    32.3379     7.7172     8.7130     <2e-16   759.6828
#> theta_2    -9.8200     0.5338   -18.3958   -10.8663    -8.7738     <2e-16   248.6404
#> theta_3     4.5745     0.8275     5.5283     2.9527     6.1963   3.23e-08    24.8825
#> theta_4    41.1707     7.4473     5.5283    26.5743    55.7671   3.23e-08    24.8823
```

# Generate predicted survival measures from a parametric survival model

Computes predicted survival analysis measures and point-wise confidence
intervals using the delta method. Meant to be used after fitting
[`ee_survival_model()`](https://r-causal.github.io/deli/reference/ee_survival_model.md)
with
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md).

## Usage

``` r
survival_predictions(
  times,
  theta,
  covariance,
  distribution,
  measure = "survival",
  alpha = 0.05,
  deriv_method = "capprox",
  dx = 1e-09
)
```

## Arguments

- times:

  Numeric vector of time points for prediction.

- theta:

  Numeric vector of estimated parameters from `ee_survival_model`.

- covariance:

  Numeric covariance matrix from `vcov(m)`.

- distribution:

  Character string matching the distribution used in
  `ee_survival_model`: `"exponential"`, `"weibull"`, or `"gompertz"`.
  Any other value is an error.

- measure:

  Character string: `"survival"`, `"risk"`, `"density"`, `"hazard"`, or
  `"cumulative_hazard"`. Default `"survival"`.

- alpha:

  Numeric significance level. Default `0.05`.

- deriv_method:

  Character string for the derivative method used to build the
  delta-method Jacobian. One of `"capprox"` (central difference),
  `"fapprox"` (forward difference), `"bapprox"` (backward difference),
  or `"exact"` (forward-mode automatic differentiation). Default
  `"capprox"`. Python Delicatessen uses exact differentiation
  internally; pass `deriv_method = "exact"` to reproduce it with exact
  derivatives and no step-size tuning. See
  [`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md).

- dx:

  Numeric step size for the finite-difference methods; ignored when
  `deriv_method = "exact"`. Default `1e-9`. Must be a single positive
  finite number, which is checked whichever `deriv_method` is in force.
  The step is absolute and is floored at the floating-point resolution
  of each estimate, so a large parameter magnitude cannot silently
  reduce it to nothing; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

## Value

A data frame with columns: `time`, `predicted`, `variance`, `lower`,
`upper`.

## Gompertz distribution

The Gompertz survival and hazard follow the
[`ee_survival_model()`](https://r-causal.github.io/deli/reference/ee_survival_model.md)
parameterization,

\$\$S(t) = \exp\left(-\frac{\lambda}{\gamma}\left(e^{\gamma t} -
1\right) \right), \qquad h(t) = \lambda e^{\gamma t}.\$\$

This is a deliberate divergence from Python Delicatessen, whose
`survival_predictions` branches only on the exponential distribution and
otherwise applies the Weibull formulas, so a `"gompertz"` request there
silently returns Weibull values.

## Examples

``` r
# Weibull survival times for 45 women with breast cancer, with no covariates.
# The default rootSolve solver does not converge here, so nleqslv is used.
psi <- function(theta) {
  ee_survival_model(
    theta,
    time = breast_cancer$times,
    event = breast_cancer$delta,
    distribution = "weibull"
  )
}
m <- m_estimate(
  stacked_equations = psi,
  init = c(0.1, 0.1),
  solver = "nleqslv"
)

# The survival function at three follow-up times, with delta-method
# confidence intervals. Observed times run from 5 to 225 months.
survival_predictions(
  times = c(50, 100, 150),
  theta = coef(m),
  covariance = vcov(m),
  distribution = "weibull",
  measure = "survival"
)
#>   time predicted    variance     lower     upper
#> 1   50 0.7209873 0.002501555 0.6229586 0.8190160
#> 2  100 0.5420930 0.004593458 0.4092564 0.6749297
#> 3  150 0.4133109 0.006038501 0.2610065 0.5656154
```

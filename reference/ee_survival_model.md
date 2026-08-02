# Estimating equation for parametric survival models

Returns a p-by-n matrix of estimating equation contributions for
parametric survival models. Supports exponential, Weibull, and Gompertz
distributions.

## Usage

``` r
ee_survival_model(theta, time, event, distribution)
```

## Arguments

- theta:

  Numeric vector of distribution parameters. For exponential, a single
  parameter (lambda). For Weibull and Gompertz, two parameters (lambda,
  gamma).

- time:

  Numeric vector of n observed (possibly censored) times.

- event:

  Numeric vector of n event indicators (1 = event, 0 = censored).

- distribution:

  Character string: `"exponential"`, `"weibull"`, or `"gompertz"`.

## Value

A p-by-n matrix where p is the number of parameters. The row is named
`lambda` for the exponential distribution, and the rows are named
`lambda` and `gamma` for the Weibull and Gompertz distributions.

## Details

The estimating equations are based on the score equations of the
corresponding parametric survival model, accounting for right censoring.
For event observations, the contribution comes from the log-density; for
censored observations, the contribution comes from the log-survival
function.

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

# The first parameter is the scale, the second the shape. A shape near 1
# means the Weibull fits about as well as the simpler exponential.
coef(m)
#>      lambda       gamma 
#> 0.009509931 0.904399702 
```

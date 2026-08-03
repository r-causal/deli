# Function-level predicted survival measures from an AFT model

Computes predicted survival analysis measures and point-wise confidence
intervals from an accelerated failure time model for a single covariate
pattern across a set of time points. The point estimates mirror
[`aft_predictions_individual()`](https://r-causal.github.io/deli/reference/aft_predictions_individual.md);
the variance is obtained with the delta method (see
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md))
and Wald-type intervals are formed on the resulting standard errors.
Meant to be used after fitting
[`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) with
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md),
typically to draw a measure and its confidence band over time.

## Usage

``` r
aft_predictions_function(
  X,
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

- X:

  Numeric 1-by-b design matrix giving a single covariate pattern. More
  than one row is an error, since each pattern has its own variance.

- times:

  Numeric vector of time points for prediction.

- theta:

  Numeric vector of estimated parameters from `ee_aft`.

- covariance:

  Numeric covariance matrix from `vcov(m)`.

- distribution:

  Character string matching the distribution used in `ee_aft`.

- measure:

  Character string: `"survival"`, `"risk"`, `"density"`, `"hazard"`, or
  `"cumulative_hazard"`. Default `"survival"`.

- alpha:

  Numeric significance level. Default `0.05` (95% CIs).

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

A data frame with one row per time point and columns: `time`,
`predicted`, `variance`, `lower`, `upper`.

## Details

The survival and hazard for a covariate pattern \\X\\ are

\$\$S(t) = S\_{\epsilon}\left( \frac{\log(t) - X \beta^T}{\sigma}
\right)\$\$ \$\$h(t) = (\sigma t)^{-1} h\_{\epsilon}\left(
\frac{\log(t) - X \beta^T}{\sigma} \right)\$\$

where \\S\_{\epsilon}\\ and \\h\_{\epsilon}\\ are the error survival and
hazard functions for the chosen `distribution`. The requested `measure`
is derived from these through
[`convert_survival_measures()`](https://r-causal.github.io/deli/reference/convert_survival_measures.md).

## Length-one times

A single time point is supported and returns a one-row data frame equal
to the corresponding row of a multi-time call. This is a deliberate
improvement over Python Delicatessen, whose `aft_predictions_function`
raises on a scalar or length-one `times` because it takes the diagonal
of a scalar delta-method covariance.

## Examples

``` r
# Weibull AFT fit, then a survival curve for one covariate pattern
set.seed(1)
n <- 200
x <- rbinom(n, 1, 0.5)
Xd <- cbind(1, x)
eps <- log(-log(runif(n)))
t_event <- exp(2 + 0.5 * x + 0.8 * eps)
t_censor <- rexp(n, rate = 0.02)
t_obs <- pmin(t_event, t_censor)
delta <- as.numeric(t_event <= t_censor)

psi <- function(theta) {
  ee_aft(theta, X = Xd, time = t_obs, event = delta, distribution = "weibull")
}
m <- m_estimate(
  stacked_equations = psi,
  init = c(mean(log(t_obs)), 0, 0),
  solver = "nleqslv"
)

aft_predictions_function(
  X = matrix(c(1, 1), nrow = 1), times = c(5, 10, 20),
  theta = coef(m), covariance = vcov(m),
  distribution = "weibull", measure = "risk"
)
#>   time predicted     variance     lower     upper
#> 1    5 0.2660215 0.0012661548 0.1962800 0.3357631
#> 2   10 0.5515778 0.0016837710 0.4711531 0.6320026
#> 3   20 0.8750468 0.0007048205 0.8230128 0.9270808
```

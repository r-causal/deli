# Predicted survival measures from an AFT model

Computes individual-level predicted survival measures from an
accelerated failure time model at specified time points. Meant to be
used after fitting
[`ee_aft()`](https://r-causal.github.io/deli/reference/ee_aft.md) with
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md).

## Usage

``` r
aft_predictions_individual(X, times, theta, distribution, measure = "survival")
```

## Arguments

- X:

  Numeric n-by-b design matrix of covariate values.

- times:

  Numeric vector of time points for prediction.

- theta:

  Numeric vector of estimated parameters from `ee_aft`.

- distribution:

  Character string matching the distribution used in `ee_aft`.

- measure:

  Character string: `"survival"`, `"risk"`, `"density"`, `"hazard"`, or
  `"cumulative_hazard"`. Default `"survival"`.

## Value

A data frame with n rows and one column per time point.

## Examples

``` r
# Weibull AFT fit, then individual-level survival for the first four people
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

# Rows are individuals and columns are the requested times. The first two
# people share a covariate pattern, as do the third and fourth, so their
# predictions agree.
aft_predictions_individual(
  X = Xd[1:4, ], times = c(5, 10, 20),
  theta = coef(m), distribution = "weibull", measure = "survival"
)
#>         t_5      t_10       t_20
#> 1 0.5988246 0.2645382 0.03179688
#> 2 0.5988246 0.2645382 0.03179688
#> 3 0.7339785 0.4484222 0.12495318
#> 4 0.7339785 0.4484222 0.12495318
```

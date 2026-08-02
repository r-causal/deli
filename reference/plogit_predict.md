# Predicted survival measures from a pooled logistic regression model

Computes predicted survival analysis measures from a pooled logistic
regression model at specified time points. Meant to be used after
fitting
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
with
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md).

## Usage

``` r
plogit_predict(
  theta,
  time,
  event,
  X,
  S = NULL,
  times_to_predict = NULL,
  measure = "survival",
  unique_times = NULL
)
```

## Arguments

- theta:

  Numeric vector of estimated parameters from `ee_plogit`.

- time:

  Numeric vector of n observed (possibly censored) times (same as in
  `ee_plogit`).

- event:

  Numeric vector of n event indicators (same as in `ee_plogit`).

- X:

  Numeric n-by-b design matrix for covariates.

- S:

  Optional time design matrix, with one row per unit-time interval up to
  the maximum observed time. Default `NULL` uses disjoint indicators.

- times_to_predict:

  Optional numeric vector of specific times to predict at. Default
  `NULL` returns all time steps.

- measure:

  Character string: `"survival"`, `"risk"`, `"density"`, `"hazard"`, or
  `"cumulative_hazard"`. Default `"survival"`.

- unique_times:

  Optional numeric vector of unique event times. Default `NULL`. When
  `S` is supplied it may only agree with the time grid the function
  builds; see the section on a supplied time design matrix.

## Value

A K-by-n matrix (or selected subset of rows if `times_to_predict` is
specified).

## A supplied time design matrix

A supplied `S` models time parametrically over the unit-time intervals
from one to the maximum observed time, one row of `S` per interval. That
grid is the function's to build, and the two arguments that describe it
have to agree with it rather than replace it: `nrow(S)` counts its steps
and `unique_times`, when supplied, names them. A mismatch in either is
an error.

Neither can be honored in place of the built grid, because the grid is
also the binning of the person-periods
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
solves on, so predictions on any other grid would come from coefficients
that were never fitted to it. A maximum observed time falling between
two whole times names no further whole interval, so the grid stops at
the last whole one and an `S` sized past it is an error as well.
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
validates both arguments the same way, so a grid the equation refuses is
not one predictions come back from.

This is a deliberate divergence from Python Delicatessen, which
documents `unique_times` as ignored when a time design matrix is
supplied. An ignored argument leaves a caller believing the grid was
theirs to choose, so deli validates it instead.

## Examples

``` r
# Bladder tumor recurrence, fit with disjoint indicators for time
W <- cbind(
  novel = collett_bladder$treat - 1,
  as.matrix(collett_bladder[, c("init", "size")])
)
k <- length(unique(collett_bladder$time[collett_bladder$delta == 1]))

psi <- function(theta) {
  ee_plogit(
    theta,
    X = W,
    time = collett_bladder$time,
    event = collett_bladder$delta
  )
}
m <- m_estimate(
  stacked_equations = psi,
  init = c(rep(0, ncol(W)), -4, rep(0, k - 1))
)

# Rows are the requested times and columns are individuals, so this shows
# disease-free survival at 12, 24, and 59 months for the first five people.
plogit_predict(
  coef(m),
  time = collett_bladder$time,
  event = collett_bladder$delta,
  X = W,
  times_to_predict = c(12, 24, 59),
  measure = "survival"
)[, 1:5]
#>           [,1]      [,2]      [,3]      [,4]       [,5]
#> [1,] 0.6310973 0.5882252 0.5534843 0.6310973 0.29096680
#> [2,] 0.5302241 0.4810942 0.4422455 0.5302241 0.18117374
#> [3,] 0.3926549 0.3402523 0.3005634 0.3926549 0.08076869
```

# Estimating equation for pooled logistic regression

Returns a p-by-n matrix of estimating equation contributions for pooled
logistic regression with discrete-time survival data. This
implementation does not require creation of a long data set.

## Usage

``` r
ee_plogit(
  theta,
  X,
  time,
  event,
  S = NULL,
  unique_times = NULL,
  weights = NULL,
  offset = NULL
)
```

## Arguments

- theta:

  Numeric vector of length `b + p_s`, where `b` is the number of
  covariate columns in `X`. When `S` is supplied, `p_s` is `ncol(S)`;
  when `S = NULL`, `p_s` is `K`, the number of unique event times.

- X:

  Numeric n-by-b design matrix for baseline covariates.

- time:

  Numeric vector of n observed (possibly censored) times.

- event:

  Numeric vector of n event indicators (1 = event, 0 = censored).

- S:

  Optional time design matrix with K rows (one per time step) and p_s
  columns. Default `NULL` uses disjoint indicators for unique event
  times. When supplied, time is modeled over the unit-time intervals
  from one to the maximum observed time, so K has to be the number of
  those intervals.

- unique_times:

  Optional numeric vector of unique event times. Default `NULL`. When
  `S = NULL` it names the time steps the disjoint indicators are built
  for. When `S` is supplied the grid is the unit-time intervals up to
  the maximum observed time, and that grid is also the binning of the
  person-periods the equation is solved on, so `unique_times` may only
  agree with it: a value equal to it is accepted and changes nothing,
  and any other value is an error rather than the silently ignored
  argument Python Delicatessen documents.
  [`plogit_predict()`](https://r-causal.github.io/deli/reference/plogit_predict.md)
  validates it the same way.

- weights:

  Optional numeric vector of n weights, or an n-by-K matrix of
  time-varying weights with one column per time interval (K must equal
  the number of unit-time intervals). Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A `(b + p_s)`-by-n matrix.

## Examples

``` r
# Bladder tumor recurrence, comparing the novel treatment to placebo while
# adjusting for the number and size of the initial tumors.
W <- cbind(
  novel = collett_bladder$treat - 1,
  as.matrix(collett_bladder[, c("init", "size")])
)

# Time is modeled with disjoint indicators, one per distinct event time,
# which ee_plogit builds by default.
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

# The first three parameters are the covariate coefficients, which
# approximate log hazard ratios. The rest describe the baseline hazard over
# time: the first of them is the log-odds of an event at the earliest event
# time, and each one after it is that time's departure from it.
summary(m, subset = 1:3)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 86
#> Parameters: 24
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1    -0.5479     0.3298    -1.6613    -1.1944     0.0985     0.0967     3.3710
#> theta_2     0.2598     0.0824     3.1535     0.0983     0.4213    0.00161     9.2757
#> theta_3     0.0735     0.0931     0.7898    -0.1089     0.2560       0.43     1.2188
```

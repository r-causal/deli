# Generate observation weights

Returns observation weights. If `weights` is `NULL`, returns a vector of
ones (equal weighting). A numeric vector is validated against `n` and
returned. A weight matrix (one row per observation, one column per time
interval, used by
[`ee_plogit()`](https://r-causal.github.io/deli/reference/ee_plogit.md)
for time-varying weights) is validated on its row count and returned
with its dimensions preserved; the column-count check is left to the
caller, which alone knows the number of intervals.

## Usage

``` r
generate_weights(n, weights = NULL)
```

## Arguments

- n:

  Integer number of observations.

- weights:

  Numeric vector of weights, an n-row weight matrix, or `NULL` for equal
  weights.

## Value

A numeric vector of length `n`, or the supplied n-row matrix.

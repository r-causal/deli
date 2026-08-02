# Check that a LASSO approximation epsilon is non-negative

Validates the `epsilon` argument of the approximate LASSO estimating
equations directly, so a negative value is rejected with a message that
names `epsilon` rather than surfacing the downstream bridge-penalty
message phrased in terms of `gamma`. Mirrors Python Delicatessen, which
validates `epsilon` up front. Zero is permitted and falls through to the
bridge penalty's non-differentiability warning.

## Usage

``` r
check_epsilon(epsilon)
```

## Arguments

- epsilon:

  The approximation parameter supplied by the caller.

## Value

Invisible `NULL`. Raises an error if `epsilon` is negative.

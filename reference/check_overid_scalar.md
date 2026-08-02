# Check that an over-identification control is a single number

Validates that a GMM over-identification control (`overid_maxiter` or
`overid_tolerance`) is a single, non-missing number. Non-positive values
are permitted: they encode degenerate but supported settings that match
Python Delicatessen, so they are not rejected here.

## Usage

``` r
check_overid_scalar(value, arg)
```

## Arguments

- value:

  The control value supplied to
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md).

- arg:

  The argument name, used in the error message.

## Value

Invisible `NULL`. Raises an error if the value is invalid.

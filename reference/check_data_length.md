# Check that a data vector matches the observation count

Validates that a data argument (`y`, `offset`, `q_eval`, `delta`, and
the like) has one value per observation, mirroring the broadcast errors
Python Delicatessen raises when an argument length does not match the
data. Tangent containers are measured on the primal, because the
response reaching a regression estimating equation can be a
`PrimalTangentArray` whose
[`length()`](https://rdrr.io/r/base/length.html) is not the observation
count.

## Usage

``` r
check_data_length(x, n, arg)
```

## Arguments

- x:

  The data argument to validate.

- n:

  The number of observations the argument must match.

- arg:

  The argument name, used in the error message.

## Value

Invisible `NULL`. Raises an error if the length does not match.

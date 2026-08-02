# Check that penalty and center have valid shapes

Validates that `penalty` and `center` are either length 1 or the same
length as `theta`, and that all penalty values are non-negative.

## Usage

``` r
check_penalty_shape(theta, penalty, center)
```

## Arguments

- theta:

  Numeric vector of parameters.

- penalty:

  Numeric penalty term (scalar or vector).

- center:

  Numeric center for penalty (scalar or vector).

## Value

Invisible `NULL`. Raises an error if shapes are invalid.

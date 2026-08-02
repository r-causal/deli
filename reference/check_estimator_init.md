# Check that init is a non-empty numeric vector

Validates that `init` is numeric and contains at least one value.

## Usage

``` r
check_estimator_init(init)
```

## Arguments

- init:

  The initial parameter vector supplied to an estimator.

## Value

Invisible `NULL`. Raises an error if `init` is invalid.

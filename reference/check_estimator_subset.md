# Check that subset holds valid parameter indices

Validates that `subset` is either `NULL` or a vector of whole-number,
1-based parameter indices within `1:n_params`.

## Usage

``` r
check_estimator_subset(subset, n_params)
```

## Arguments

- subset:

  The parameter subset supplied to an estimator.

- n_params:

  The number of parameters, taken from `length(init)`.

## Value

Invisible `NULL`. Raises an error if `subset` is invalid.

# Check that survival data is valid

Validates that event indicators are 0 or 1 (ignoring NAs) and that
observation times are positive (ignoring NAs).

## Usage

``` r
check_survival_data_valid(delta, time)
```

## Arguments

- delta:

  Numeric vector of event indicators (0 or 1).

- time:

  Numeric vector of observation times.

## Value

Invisible `NULL`. Raises an error if data is invalid.

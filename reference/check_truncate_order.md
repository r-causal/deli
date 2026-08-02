# Check that truncation bounds are in ascending order

Validates that a length-2 `truncate` vector has its lower bound no
greater than its upper bound, mirroring the check Python Delicatessen
performs before clipping propensity scores.

## Usage

``` r
check_truncate_order(truncate)
```

## Arguments

- truncate:

  Length-2 numeric vector `c(lower, upper)`.

## Value

Invisible `NULL`. Raises an error if bounds are out of order.

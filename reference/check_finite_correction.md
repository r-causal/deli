# Check that finite_correction names a supported correction

Validates that `finite_correction` is either `NULL` or a supported
correction string. `"HC1"` is the only supported non-`NULL` value,
matching
[`finite_sample_correction()`](https://r-causal.github.io/deli/reference/finite_sample_correction.md).

## Usage

``` r
check_finite_correction(finite_correction)
```

## Arguments

- finite_correction:

  The finite-sample correction supplied to an estimator.

## Value

Invisible `NULL`. Raises an error if the value is unsupported.

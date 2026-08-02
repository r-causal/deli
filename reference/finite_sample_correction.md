# Apply finite-sample correction to the meat matrix

Applies the HC1 correction: \\meat \times n / (n - p)\\.

## Usage

``` r
finite_sample_correction(meat, n, p, adjustment = NULL)
```

## Arguments

- meat:

  A p-by-p meat matrix.

- n:

  Integer number of observations.

- p:

  Integer number of parameters.

- adjustment:

  Character string or `NULL`. Currently only `"HC1"` is supported.

## Value

The corrected meat matrix.

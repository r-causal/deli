# Compute the meat matrix

Computes the meat matrix as the cross-product of the estimating equation
evaluations: \\EE \times EE^T\\.

## Usage

``` r
compute_meat(evaluations)
```

## Arguments

- evaluations:

  A p-by-n matrix of estimating equation evaluations, where p is the
  number of parameters and n is the number of observations.

## Value

A p-by-p meat matrix.

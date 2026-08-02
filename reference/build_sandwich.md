# Build the sandwich variance estimator

Combines bread and meat matrices into the sandwich: \\B^{-1} M
(B^{-1})^T\\.

## Usage

``` r
build_sandwich(bread, meat, allow_pinv = TRUE, call = rlang::caller_env())
```

## Arguments

- bread:

  A bread matrix with one row per estimating equation and one column per
  parameter. It is p-by-p for an M-estimation system, which has one
  equation per parameter, and n_eqs-by-p for an over-identified GMM
  system, whose rectangular bread has no inverse and is pseudo-inverted
  instead.

- meat:

  An n_eqs-by-n_eqs meat matrix, square in the estimating equations
  whether or not the bread is.

- allow_pinv:

  Logical. If `TRUE` (default), uses the pseudo-inverse when the bread
  matrix cannot be inverted. When `FALSE`, a bread that has no inverse
  raises an error carrying the class `deli_bread_not_invertible`.

- call:

  The frame to report that error against.

## Value

A p-by-p sandwich covariance matrix, or `NULL` if the bread contains
`NA` values.

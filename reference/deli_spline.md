# Generate polynomial spline basis terms

Generates polynomial spline terms for a numeric vector at pre-specified
knot locations. Default is restricted (natural) cubic splines, but
unrestricted splines with different polynomial terms can also be
generated.

This function mirrors
[`spline()`](https://rdrr.io/r/stats/splinefun.html) in Python
delicatessen, so code translated from Python can keep its shape. Despite
the name it has nothing to do with
[`stats::spline()`](https://rdrr.io/r/stats/splinefun.html), which
interpolates a curve through data points and returns the interpolated
values. `deli_spline()` builds a truncated power basis, a matrix of
design columns to be used as covariates in a regression. The `deli_`
prefix is there so that the two names do not collide.

No base R function returns these truncated power columns, but the
`splines` package, which installs with every copy of R, builds the same
spline spaces in a B-spline parameterization.
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) is the counterpart
to `restricted = FALSE`: for knots `k`,
`cbind(1, x, x^2, x^3, deli_spline(x, k, restricted = FALSE))` and
`cbind(1, splines::bs(x, knots = k))` span the same space and give
identical fitted values.
[`splines::ns()`](https://rdrr.io/r/splines/ns.html) is the nearest
counterpart to `restricted = TRUE` but not an exact one, since it also
constrains the fit to be linear beyond the boundary knots and so spans a
subspace of the columns here. `deli_spline()` is what deli offers for
parity with Python delicatessen, for the truncated power form itself,
and because
[`additive_design_matrix()`](https://r-causal.github.io/deli/reference/additive_design_matrix.md)
builds on it.

## Usage

``` r
deli_spline(x, knots, power = 3, restricted = TRUE, normalized = FALSE)
```

## Arguments

- x:

  Numeric vector of observed values.

- knots:

  Numeric vector of knot locations. Should be between the min and max of
  `x`.

- power:

  Numeric power for the spline terms. Default `3` (cubic).

- restricted:

  Logical. If `TRUE` (default), generate restricted (natural) splines.
  If `FALSE`, generate unrestricted splines.

- normalized:

  Logical. If `TRUE`, divide spline terms by the range of knots (largest
  minus smallest) raised to `power`. With a single knot the range is
  zero, so the divisor is that knot raised to `power` instead. Default
  `FALSE`.

## Value

A matrix with `length(x)` rows. Number of columns is `length(knots)` for
unrestricted or `length(knots) - 1` for restricted.

## Details

Unrestricted splines for knot \\k\\ are: \$\$s_k(X) = I(X \> k) (X -
k)^a\$\$

Restricted (natural) splines subtract the last spline term: \$\$r_k(X) =
s_k(X) - s_K(X)\$\$ where \\K\\ is the largest knot. Restricted splines
return one fewer column than the number of knots.

## Examples

``` r
# Restricted quadratic splines at four knots, so three basis columns
s <- deli_spline(1:59, knots = c(10, 20, 30, 40), power = 2)
dim(s)
#> [1] 59  3

# Each term stays at zero below its knot and grows above it
s[c(5, 15, 25, 35, 45, 55), ]
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]   25    0    0
#> [3,]  225   25    0
#> [4,]  625  225   25
#> [5,] 1200  600  200
#> [6,] 1800 1000  400
```

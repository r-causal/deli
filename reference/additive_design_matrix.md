# Build an additive design matrix for GAMs

Constructs the expanded design matrix for generalized additive models by
appending spline basis terms according to per-column specifications.
Each column in `X` keeps its linear term; columns with a non-`NULL`
specification get additional spline basis columns appended.

This function mirrors `additive_design_matrix()` in Python delicatessen,
so code translated from Python can keep its shape. There is no base R
equivalent of this construction, so this is the interface for it in deli
as well.

## Usage

``` r
additive_design_matrix(X, specifications, return_penalty = FALSE)
```

## Arguments

- X:

  Numeric matrix (n-by-b) of input covariates.

- specifications:

  A list of length `b` (number of columns in `X`). Each element is
  either `NULL` (no spline for that column) or a list with:

  knots

  :   Numeric vector of knot locations (required).

  natural

  :   Logical, generate restricted splines? Default `TRUE`.

  power

  :   Numeric power for spline. Default `3`.

  penalty

  :   Numeric penalty for spline terms. Default `0`.

  normalized

  :   Logical, normalize spline terms? Default `FALSE`.

- return_penalty:

  Logical. If `TRUE`, return a list with both the design matrix and the
  penalty vector. Default `FALSE`.

## Value

If `return_penalty = FALSE`, a numeric matrix. If `TRUE`, a list with
elements `X` (the design matrix) and `penalty` (numeric vector).

## Examples

``` r
set.seed(42)
X <- cbind(rnorm(50), rnorm(50))

# The first column stays linear; the second also gets penalized splines
specs <- list(NULL, list(knots = c(-1, 0, 1), penalty = 5))
out <- additive_design_matrix(X, specs, return_penalty = TRUE)

# Linear terms are unpenalized, spline terms carry the requested penalty
out$penalty
#> [1] 0 0 5 5

dim(out$X)
#> [1] 50  4
```

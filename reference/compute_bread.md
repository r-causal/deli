# Compute the bread matrix

Computes the bread matrix for the empirical sandwich variance estimator.
The bread is the negative Jacobian of the summed estimating equations.
It is returned unscaled: the callers that assemble a sandwich divide it
by the number of observations, and the meat with it.

## Usage

``` r
compute_bread(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-09,
  summed_equations = NULL
)
```

## Arguments

- stacked_equations:

  A function that takes a numeric vector `theta` and returns a p-by-n
  matrix of estimating equation contributions.

- theta:

  Numeric vector of parameter estimates.

- deriv_method:

  Character string for the derivative method. One of `"capprox"`
  (central), `"fapprox"` (forward), or `"bapprox"` (backward).

- dx:

  Numeric step size (default `1e-9`). The step is absolute, floored at
  the floating-point resolution of each estimate; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

- summed_equations:

  A function of `theta` returning the length-p vector of row sums of
  `stacked_equations` at `theta`, or `NULL` (default) to derive that
  reduction from the full p-by-n return. See
  [`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md)
  for what a supplied reduction saves and what it must satisfy.

  Anything that is neither `NULL` nor a function, and any return at
  `theta` that is not numeric or holds fewer values than there are
  parameters, raises an error carrying the class
  `deli_summed_equations_error`. The values themselves are taken on
  trust: this function evaluates the estimating equations nowhere, so it
  has nothing to compare them against, and a reduction that sums some
  other system returns the Jacobian of that other system with nothing to
  say so.

## Value

The negated Jacobian of the summed estimating equations, with one row
per estimating equation and one column per parameter. That is p-by-p for
an M-estimation system, which has one equation per parameter, and
n_eqs-by-p for an over-identified GMM system, whose rectangular bread
[`build_sandwich()`](https://r-causal.github.io/deli/reference/build_sandwich.md)
pseudo-inverts. No scaling is applied here; the division by n that puts
the bread on the mean scale belongs to the callers that assemble a
sandwich,
[`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md)
and
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md).

A bread holding `NA` is returned as it stands, alongside a warning
carrying the class `deli_bread_na`. What to do about it is the caller's,
and the two callers differ: a fit records no variance and says so, while
[`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md)
has nothing but the matrix to return and fails.

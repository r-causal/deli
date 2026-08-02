# Compute the empirical sandwich variance estimator

Computes the empirical sandwich variance estimator directly from a set
of estimating equations and a vector of parameter estimates. Unlike
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md),
this function does not solve for the parameters; it assumes that `theta`
is already the root of the estimating equations and only assembles the
sandwich covariance at that point.

## Usage

``` r
compute_sandwich(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-09,
  allow_pinv = TRUE,
  finite_correction = NULL,
  summed_equations = NULL,
  check_summed_equations = TRUE
)
```

## Arguments

- stacked_equations:

  A function that takes a numeric vector `theta` and returns a p-by-n
  matrix of estimating equation contributions, where p is the number of
  parameters and n is the number of observations. A list of one element
  per equation, each holding that equation's contributions across the
  observations, is accepted as well.
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  has no support for that shape, so an estimating function written in it
  reaches a variance through this entry point and a Jacobian through
  [`compute_bread()`](https://r-causal.github.io/deli/reference/compute_bread.md).
  The list form holds exactly one element per parameter, each of one
  length, so an over-identified system reaches this function as a matrix
  rather than as a list.

- theta:

  Numeric vector of parameter estimates. This function assumes `theta`
  is the root of `stacked_equations`; it does not solve for it.

- deriv_method:

  Character string selecting the method used to build the bread
  Jacobian. One of `"capprox"` (central difference), `"fapprox"`
  (forward difference), `"bapprox"` (backward difference), or `"exact"`
  (forward-mode automatic differentiation). Default `"capprox"`.

  This default differs from Python delicatessen, whose
  `compute_sandwich` defaults to `"approx"`, a forward difference
  computed through SciPy's `approx_fprime`. deli does not replicate the
  SciPy `"approx"` path; its `"fapprox"` is the hand-implemented forward
  difference. Code ported across the two libraries should set
  `deriv_method` explicitly rather than rely on the default.

- dx:

  Numeric step size for the finite-difference methods; ignored when
  `deriv_method = "exact"`. A small value is recommended, since large
  steps can give poor approximations. Default `1e-9`. Must be a single
  positive finite number, which is checked whichever `deriv_method` is
  in force. The step is absolute and is floored at the floating-point
  resolution of each estimate, so a large parameter magnitude cannot
  silently reduce it to nothing; see
  [`approx_differentiation()`](https://r-causal.github.io/deli/reference/approx_differentiation.md).

- allow_pinv:

  Logical. When `TRUE` (default), the Moore-Penrose pseudo-inverse is
  used where the bread matrix cannot be solved, which is the case for
  the rectangular bread of an over-identified system as well as for a
  square one whose solve fails; when `FALSE`, a bread with no inverse
  raises an error carrying the class `deli_bread_not_invertible`.

  What the two settings choose between is what to do when the solve
  fails, and nothing else: a bread the solve inverts is inverted under
  either. An ill-conditioned bread that base R returns an inverse for is
  one of those, however far its condition number runs, so
  `allow_pinv = FALSE` is not a conditioning test and a caller who wants
  one applies it to the returned matrix.

- finite_correction:

  Character string or `NULL`. Finite-sample correction applied to the
  meat matrix. `NULL` (default) applies no correction; `"HC1"` rescales
  the meat by \\n / (n - p)\\, where p is the number of parameters.

- summed_equations:

  A function that takes a numeric vector `theta` and returns the
  length-p vector of row sums of `stacked_equations` at `theta`, or
  `NULL` (default) to derive those sums from the full p-by-n return.

  The bread is the Jacobian of the summed estimating equations, so each
  of the one or two perturbed evaluations it makes per parameter is
  reduced to one value per equation as soon as it is built. Deriving the
  reduction builds the whole p-by-n matrix 2p times for arithmetic that
  is linear in it; an estimating function whose sums have a closed form,
  such as the \\X^T r\\ of a regression score, can supply them and hand
  the bread only what it uses. The meat is unaffected either way, since
  it needs the per-observation contributions and takes the one full
  evaluation it always took.

  Under `deriv_method = "exact"` the reduction is called with a
  tangent-carrying `theta`, so it must be written in operations that
  carry derivatives: `t(X) %*% r` does, and
  [`base::crossprod()`](https://rdrr.io/r/base/crossprod.html) does not.
  See
  [`auto_differentiation()`](https://r-causal.github.io/deli/reference/auto_differentiation.md)
  for which operations carry a tangent and where.

  An argument that is neither `NULL` nor a function raises an error
  carrying the class `deli_summed_equations_error`, whatever
  `check_summed_equations` says. So does a return at `theta` that is not
  numeric, or one holding fewer values than there are parameters, both
  of which
  [`compute_bread()`](https://r-causal.github.io/deli/reference/compute_bread.md)
  reads before it differentiates anything. That the return holds exactly
  one value per estimating equation is read by the comparison below,
  since only a call that has evaluated the estimating functions knows
  how many of them there are.

- check_summed_equations:

  Logical. When `TRUE` (default) and `summed_equations` was supplied,
  its value at `theta` is compared against the row sums of the one full
  evaluation the meat is built from, and a disagreement raises an error
  carrying the class `deli_summed_equations_disagree`. The comparison
  costs one call to `summed_equations` and one reduction of a matrix the
  call already holds. Set it to `FALSE` to skip the comparison, which
  leaves a reduction that sums some other system to return a matrix with
  the shape of a covariance and no claim to be one.

  The comparison is made at `theta`, which the caller states is the
  root, so both quantities are at rounding there. That is enough to
  catch a reduction of some other system and not enough to catch a
  reduction that is a multiple of the right one, which agrees at a root
  and yields that multiple of the right bread.

## Value

A p-by-p covariance matrix on the asymptotic scale. The bread and meat
are each divided by n internally, so the returned matrix is the variance
that corresponds to the standard deviation. Dividing it by the number of
observations gives the standard-error-scale variance, whose
square-rooted diagonal is the vector of standard errors.

A covariance matrix is the only thing this function returns. Where the
bread has no inverse, and so no sandwich can be assembled from it, the
call raises an error carrying the class `deli_bread_not_invertible`
rather than returning something that has to be tested for. That covers a
bread holding `NA`, and, under `allow_pinv = FALSE`, a rectangular bread
and one the solve could not invert.
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
makes the other choice from the same matrices: a fit whose bread holds
`NA` warns and comes back with no variance, since it still carries the
estimates.

## Details

The sandwich is built from a bread matrix and a meat matrix. The bread
is the negative Jacobian of the summed estimating equations, and the
meat is the cross-product of the equation evaluations. Each is scaled by
`1/n` internally, so the returned matrix is on the asymptotic scale (see
`Value`). The bread Jacobian is obtained either by finite differences or
by forward-mode automatic differentiation.

## References

Boos DD, & Stefanski LA. (2013). M-estimation (estimating equations). In
Essential Statistical Inference (pp. 297-337). Springer, New York, NY.

## See also

[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md)
and
[`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md),
which solve for `theta` and report this variance internally, and
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md)
for the variance of a transformation of the parameters.

## Examples

``` r
# A generic data set for estimating a mean and variance
y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)

# The mean and variance are the roots of ee_mean_variance, so they can be
# computed directly rather than solved for
theta <- c(mean(y), stats::var(y) * (length(y) - 1) / length(y))

# Wrap the built-in estimating equation as a function of theta alone
psi <- function(theta) ee_mean_variance(theta, y = y)

# compute_sandwich() returns the asymptotic-scale variance, so dividing by
# n puts it on the standard-error scale
sandwich <- compute_sandwich(psi, theta = theta) / length(y)
sandwich
#>           [,1]      [,2]
#> [1,] 0.1975308 0.2057613
#> [2,] 0.2057613 0.4883401

# The diagonal square roots are the standard errors
sqrt(diag(sandwich))
#> [1] 0.4444444 0.6988134

# The bread only ever needs the sums of the estimating equations, so an
# equation whose sums have a closed form can supply them directly
summed <- function(theta) {
  c(
    sum(y) - length(y) * theta[1],
    sum((y - theta[1])^2) - length(y) * theta[2]
  )
}
compute_sandwich(psi, theta = theta, summed_equations = summed) / length(y)
#>           [,1]      [,2]
#> [1,] 0.1975308 0.2057613
#> [2,] 0.2057613 0.4883401
```

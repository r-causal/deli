# M-Estimator

S7 class for M-estimation via solving estimating equations with
empirical sandwich variance estimation.

## Usage

``` r
MEstimator(
  stacked_equations,
  init,
  subset = NULL,
  finite_correction = NULL,
  summed_equations = NULL,
  check_summed_equations = TRUE
)
```

## Arguments

- stacked_equations:

  A function that takes a numeric vector `theta` and returns a p-by-n
  matrix of estimating equation contributions, where p is the number of
  parameters and n is the number of observations. Row names on that
  matrix name the parameters when `init` has none and every parameter is
  labeled; see
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md).

- init:

  Numeric vector of initial parameter values for the root-finding
  algorithm. Names on it label the parameters and take precedence over
  the row names of `stacked_equations`.

- subset:

  Integer vector of parameter indices to solve for, or `NULL` (default)
  to solve for all parameters. Indices are 1-based; parameters not
  listed are held fixed at their `init` values while the rest are
  solved. The equations outside the subset are set aside along with the
  parameters they estimate, so the subset parameters are the root of the
  subset equations alone and the rest of the stack has no say in where
  they land: give a three-equation linear regression stack `subset = 1L`
  and the intercept comes back as the mean of the response less what the
  slopes held at their `init` values account for, because the first
  equation on its own is the estimating equation for a mean. Held at
  zero, which is what an unset `init` usually means, they account for
  nothing and the intercept is the mean of the response itself.
  [`GMMEstimator()`](https://r-causal.github.io/deli/reference/GMMEstimator.md)
  and
  [`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
  read the argument differently, since the GMM objective sums every
  equation whether the subset lists it or not, so the same stack and the
  same `subset` give the two different values. The variance estimator
  ignores `subset`.

- finite_correction:

  Character string for finite-sample correction (e.g., `"HC1"`), or
  `NULL` (default) for no correction.

- summed_equations:

  A function that takes a numeric vector `theta` and returns the
  length-p vector of row sums of `stacked_equations` at `theta`, or
  `NULL` (default) to derive those sums from the full p-by-n return.

  A fit reduces the estimating functions to those sums everywhere but
  the meat. The solver is given the summed equations to find a root of,
  and the bread is their Jacobian, so a fit builds the whole p-by-n
  matrix once per solver evaluation and once or twice more per
  parameter, all of it for arithmetic that is linear in it. An
  estimating function whose sums have a closed form, such as the \\X^T
  r\\ of a regression score, can supply them here and give both steps
  only what they use. The meat is unaffected either way: it needs the
  per-observation contributions and takes the one full evaluation it
  always took, as does the validation at the starting values.

  Under `deriv_method = "exact"` the reduction is called with a
  tangent-carrying `theta`, so it must be written in operations that
  carry derivatives: `t(X) %*% r` does, and
  [`base::crossprod()`](https://rdrr.io/r/base/crossprod.html) does not.
  See
  [`auto_differentiation()`](https://r-causal.github.io/deli/reference/auto_differentiation.md)
  for which operations carry a tangent and where.

  Anything that is neither `NULL` nor a function is refused here, with
  an error carrying the class `deli_summed_equations_error`. So is a
  return at the estimated values that is not numeric or does not hold
  one value per estimating equation, which
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  reads.

  [`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md)
  and
  [`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md)
  do not offer it. Both build `stacked_equations` themselves from a
  formula and the equation named in `.ee`, so the estimating functions a
  reduction would have to match are ones the caller never writes; this
  is for a system assembled by hand.

- check_summed_equations:

  Logical. When `TRUE` (default) and `summed_equations` was supplied,
  its value at the estimated values is compared against the row sums of
  the one full evaluation the meat is built from, and a disagreement
  raises an error carrying the class `deli_summed_equations_disagree`.
  The comparison costs one call to `summed_equations` and one reduction
  of a matrix the fit already holds.

  What it is for is a reduction that sums some other system. The solver
  drives whatever it is given to zero, so such a reduction sends the fit
  to that other system's root and leaves the bread the Jacobian of one
  system and the meat the cross-product of another. Set it to `FALSE` to
  skip the comparison, which leaves a fit that reports estimates nobody
  solved for and a covariance with the shape of one and no claim to be
  one.

  A reduction that is a multiple of the right one is not what the
  comparison sees. It vanishes where the right one does, so the fit
  lands at the same estimates and both quantities are at rounding where
  they are compared; the bread comes back as that multiple of the right
  bread. See
  [`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md),
  whose argument of the same name reads the same comparison at a point
  the caller supplies.

## Value

An `MEstimator` S7 object. Call
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md) to
solve the estimating equations and compute the sandwich variance.

## Examples

``` r
# Estimating equations for the mean
y <- c(1, 2, 3, 4, 5)
psi <- function(theta) {
  matrix(y - theta[1], nrow = 1)
}
m <- MEstimator(stacked_equations = psi, init = 0) |>
  estimate()
coef(m)
#> theta_1 
#>       3 
summary(m)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 5
#> Parameters: 1
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1     3.0000     0.6325     4.7434     1.7604     4.2396    2.1e-06    18.8602

# The solver and the bread only ever need the sums of the estimating
# equations, so an equation whose sums have a closed form can supply them
# and leave the whole matrix to the meat
summed <- function(theta) sum(y) - length(y) * theta[1]
MEstimator(stacked_equations = psi, init = 0, summed_equations = summed) |>
  estimate() |>
  summary()
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 5
#> Parameters: 1
#> 
#>           Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> theta_1     3.0000     0.6325     4.7434     1.7604     4.2396    2.1e-06    18.8602
```

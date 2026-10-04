# Influence functions for M-Estimator

Computes the influence function values for individual observations:
\$\$IF(O_i; \theta) = B_n(\theta)^{-1} \psi(O_i; \theta)\$\$

This function mirrors `m.influence_functions()` in Python delicatessen,
so code translated from Python can keep its shape. There is no base R
accessor for influence function values, so this is the interface for
them in deli as well.

## Usage

``` r
influence_functions(object, allow_pinv = TRUE, ...)
```

## Arguments

- object:

  A fitted `MEstimator` object (after calling
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)).

- allow_pinv:

  Logical. Use pseudo-inverse if bread is singular? Default `TRUE`.

- ...:

  Not used. Must be empty, so a name that is not one of the documented
  arguments is an error rather than silently ignored.

## Value

An n-by-p matrix of influence function values, where n is the number of
observations and p is the number of parameters. Columns are named for
the parameters, as in
[`coef()`](https://r-causal.github.io/deli/reference/deli-generics.md).
The rows take whatever labels the estimating function put on the columns
of its own return, so a fit whose estimating function collapses its
contributions with
[`aggregate_efuncs()`](https://r-causal.github.io/deli/reference/aggregate_efuncs.md)
has one row per group, labeled with the group value, and every row is
that group's influence rather than an observation's.

## Examples

``` r
fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
                  model = "linear")

# One row per observation, showing its contribution to each estimate
head(influence_functions(fit))
#>       (Intercept)           wt           hp
#> [1,]  -7.88511089  1.236467576  0.009099131
#> [2,]  -3.44007774 -0.007985771  0.012831874
#> [3,] -10.16026284  1.851173715  0.011785054
#> [4,]   0.13254570  0.088097892 -0.001915625
#> [5,]   0.08947257 -0.030823617  0.002607096
#> [6,]  -0.28620413 -2.869668817  0.048707825
```

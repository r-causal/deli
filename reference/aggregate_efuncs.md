# Aggregate estimating function contributions by group

Collapses unit-level estimating function contributions into group-level
contributions for clustered or grouped data. Uses an independent working
correlation structure (summing within groups).

This function mirrors `aggregate_efuncs()` in Python delicatessen, so
code translated from Python can keep its shape. There is no base R
equivalent that operates on estimating function contributions, so this
is the interface for them in deli as well.

## Usage

``` r
aggregate_efuncs(est_funcs, group)
```

## Arguments

- est_funcs:

  A p-by-n matrix of estimating function contributions, where p is the
  number of parameters and n is the number of observations. A length-n
  vector is treated as a single parameter observed across n
  observations, matching a 1-by-n matrix.

- group:

  A vector of length n identifying the group (cluster) for each
  observation.

## Value

A p-by-m matrix, where m is the number of unique groups. Row names are
those of `est_funcs`, since the rows are the same parameters. Columns
are ordered by the sorted unique values of `group` and are labeled with
those values, as character. A factor `group` is coerced with
[`as.vector()`](https://rdrr.io/r/base/vector.html) to its character
labels before sorting, so its columns sort lexically by label rather
than by factor-level order and carry those labels; a level with no
observations contributes no column and so no label.

## Details

This function should be called inside the `psi` function after computing
unit-level estimating equations but before returning them to
[`MEstimator()`](https://r-causal.github.io/deli/reference/MEstimator.md).
This changes the effective sample size used by the empirical sandwich
variance estimator.

## Examples

``` r
# Fifty clusters of four observations, sharing a cluster-level shift in y
set.seed(42)
n <- 200
group <- rep(1:50, each = 4)
cluster_effect <- rnorm(50, sd = 2)
y <- cluster_effect[group] + rnorm(n)

psi <- function(theta) aggregate_efuncs(ee_mean(theta, y = y), group = group)

m <- m_estimate(stacked_equations = psi, init = mean(y))

# Cluster-robust standard error, larger than the naive independence version
sqrt(diag(vcov(m)))
#>   theta_1 
#> 0.3253477 
```

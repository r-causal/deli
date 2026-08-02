# Estimating equation for the percentile

Returns a 1-by-n matrix for the q-th percentile: \\\psi_i(\theta) = q -
I(Y_i \le \theta)\\.

## Usage

``` r
ee_percentile(theta, y, q)
```

## Arguments

- theta:

  Numeric vector of length 1.

- y:

  Numeric vector of observed values.

- q:

  Numeric percentile, must be in `(0, 1)`.

## Value

A 1-by-n matrix.

## Details

The derivative of this estimating equation is not defined at
\\\hat{\theta}\\, so the bread matrix and sandwich variance cannot be
used to estimate the variance. The function warns for this reason. A
direct call warns every time; a call from within
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md),
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md),
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md) or
[`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md),
each of which evaluates the estimating function many times, delivers the
warning once for the operation. It is offered for completeness but is
not generally recommended for applications.

Pass the sample quantile as `init`. The estimating function is a step
function of `theta`, so its derivative is zero wherever it is defined
and a root finder has no direction in which to search: it returns the
starting values it was given. Starting from zero therefore returns zero
rather than the quantile.

## Examples

``` r
y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
psi <- function(theta) ee_percentile(theta, y = y, q = 0.5)

# The root finder cannot move away from its starting values here, so start at
# the sample quantile, which is the solution. The fit warns once that the
# estimating equation is not differentiable, so the sandwich variance should
# not be trusted.
m <- m_estimate(stacked_equations = psi, init = median(y))
#> Warning: The estimating equation is not differentiable at `theta`. Therefore, the bread
#> matrix is not defined for finite samples, and the sandwich should not be used
#> to estimate the variance.
coef(m)
#> theta_1 
#>       3 
```

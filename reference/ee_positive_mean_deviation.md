# Estimating equations for the positive mean deviation

Returns a 2-by-n matrix for the positive mean deviation and median.

## Usage

``` r
ee_positive_mean_deviation(theta, y)
```

## Arguments

- theta:

  Numeric vector of length 2. `theta[1]` is the positive mean deviation,
  `theta[2]` is the median.

- y:

  Numeric vector of observed values.

## Value

A 2-by-n matrix, with rows named `positive_mean_deviation` and `median`.

## Details

The derivative of the estimating equation for the median is not defined
at \\\hat{\theta}\\, so the bread matrix and sandwich variance cannot be
used to estimate the variance. The function warns for this reason. A
direct call warns every time; a call from within
[`m_estimate()`](https://r-causal.github.io/deli/reference/m_estimate.md),
[`gmm_estimate()`](https://r-causal.github.io/deli/reference/gmm_estimate.md),
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md) or
[`compute_sandwich()`](https://r-causal.github.io/deli/reference/compute_sandwich.md),
each of which evaluates the estimating function many times, delivers the
warning once for the operation. It is offered for completeness but is
not generally recommended for applications.

Start `theta[2]` at the sample median and `theta[1]` at the matching
positive mean deviation, `mean(2 * (y - median(y)) * (y > median(y)))`,
and use `solver = "lm"`. The median equation is a step function of
`theta[2]`, so its derivative is zero wherever it is defined, the
finite-difference approximation the solvers work from is zero along that
row as well, and the Jacobian is singular. Neither solver recovers the
median from that, and they fail differently. The default `rootSolve`
solver cannot move at all and returns `init` unchanged for both
parameters. Whether it warns that it did not converge depends on the
starting values, so a fit can come back holding the values it was given
with nothing reported. The Levenberg-Marquardt solver drives the first
equation to zero. That equation is itself a function of `theta[2]`, so
the solver can reduce its residual by trading `theta[2]` against
`theta[1]`. Whether it makes that trade depends on the starting values:
from some it leaves the median at the value it was given, and from
others it moves the median to whatever value the trade leaves it at. It
reports no failure of its own either way, so the second outcome is
caught by the fit rather than by the solver: a median equation that the
solve left further from zero than it was at the starting values, in a
stack whose bread carries its row as zeros, warns that the estimating
equations are not solved at the returned values. Neither outcome
recovers the sample median on its own, and when the median does move,
the positive mean deviation returned belongs to that median rather than
to the sample median.

Those starting values come from measurement: forty exponential samples
per configuration, each fit with the Levenberg-Marquardt solver. Started
at the sample median with `theta[1]` at the matching deviation, the fit
returned the sample median for all forty at every one of the sizes 9,
10, 25, 40, 41, 100, and 101. Started at the sample median with
`theta[1]` at zero, it returned the sample median for about half of them
at the even sizes `n = 10` and `n = 40`, and for all forty at `n = 9`,
`n = 25`, and `n = 41`. Started with `theta[1]` at zero and `theta[2]`
half a unit, one unit, or three units to either side of the sample
median, it returned the sample median for none of the forty, at each of
those six offsets and each of the sizes 10, 25, 40, and 41.

## Examples

``` r
y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
psi <- function(theta) ee_positive_mean_deviation(theta, y = y)

# Start the median at the sample median and the deviation at the value that
# matches it, since the solver cannot search for the median, and use the
# Levenberg-Marquardt solver, which holds there. The fit warns once that the
# median estimating equation is not differentiable, so the sandwich variance
# should not be trusted.
init <- c(mean(2 * (y - median(y)) * (y > median(y))), median(y))
m <- m_estimate(stacked_equations = psi, init = init, solver = "lm")
#> Warning: The estimating equation for the median is not differentiable. Therefore, the
#> bread matrix is not defined for finite samples, and the sandwich should not be
#> used to estimate the variance.
coef(m)
#> positive_mean_deviation                  median 
#>                       2                       3 
```

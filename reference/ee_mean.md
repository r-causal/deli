# Estimating equation for the mean

Returns a 1-by-n matrix of estimating equation contributions for the
mean: \\\psi_i(\theta) = w_i (Y_i - \theta)\\.

## Usage

``` r
ee_mean(theta, y, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length 1.

- y:

  Numeric vector of observed values.

- weights:

  Optional numeric vector of weights (same length as `y`). Default
  `NULL` assigns weight 1 to all observations.

## Value

A 1-by-n matrix.

## Examples

``` r
y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
psi <- function(theta) ee_mean(theta, y = y)
m <- m_estimate(stacked_equations = psi, init = 0)
coef(m)
#> theta_1 
#>     3.4 
```

# Estimating equation for the geometric mean

Returns a 1-by-n matrix of estimating equation contributions for the
geometric mean. When `log_theta = TRUE` (default), solves \\\log(Y_i) -
\log(\theta)\\. When `log_theta = FALSE`, solves \\\log(Y_i) - \theta\\
where \\\theta\\ is the log of the geometric mean.

## Usage

``` r
ee_mean_geometric(theta, y, weights = NULL, log_theta = TRUE)
```

## Arguments

- theta:

  Numeric vector of length 1.

- y:

  Numeric vector of positive observed values.

- weights:

  Optional numeric vector of weights. Default `NULL`.

- log_theta:

  Logical. If `TRUE` (default), internally log-transforms theta. Default
  `TRUE`.

## Value

A 1-by-n matrix.

## Examples

``` r
y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
psi <- function(theta) ee_mean_geometric(theta, y = y)
m <- m_estimate(stacked_equations = psi, init = 1)
coef(m)
#>  theta_1 
#> 2.805809 
```

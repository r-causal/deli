# Estimating equations for the mean and variance

Returns a 2-by-n matrix of estimating equation contributions for the
mean and variance: \$\$\psi_i(\theta) = \begin{pmatrix} Y_i - \theta_1
\\ (Y_i - \theta_1)^2 - \theta_2 \end{pmatrix}\$\$

## Usage

``` r
ee_mean_variance(theta, y)
```

## Arguments

- theta:

  Numeric vector of length 2. `theta[1]` is the mean, `theta[2]` is the
  variance.

- y:

  Numeric vector of observed values.

## Value

A 2-by-n matrix, with rows named `mean` and `variance`.

## Examples

``` r
y <- c(1, 2, 3, 1, 4, 5, 3, 2, 6, 7)
psi <- function(theta) ee_mean_variance(theta, y = y)
m <- m_estimate(stacked_equations = psi, init = c(0, 1))
coef(m)
#>     mean variance 
#>     3.40     3.84 
```

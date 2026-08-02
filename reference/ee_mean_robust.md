# Estimating equation for the robust mean

Returns a 1-by-n matrix of estimating equation contributions for the
robust mean using the specified loss function: \\\psi_i(\theta) =
f_k(Y_i - \theta)\\.

## Usage

``` r
ee_mean_robust(theta, y, k, loss = "huber")
```

## Arguments

- theta:

  Numeric vector of length 1.

- y:

  Numeric vector of observed values.

- k:

  Numeric tuning parameter for the loss function.

- loss:

  Character string specifying the loss function. Default `"huber"`. See
  [`robust_loss_functions()`](https://r-causal.github.io/deli/reference/robust_loss_functions.md)
  for options.

## Value

A 1-by-n matrix.

## Examples

``` r
# Forty-nine standard normal observations plus one gross outlier.
set.seed(1)
y <- c(rnorm(49), 50)
psi <- function(theta) ee_mean_robust(theta, y = y, k = 1.345, loss = "huber")
m <- m_estimate(stacked_equations = psi, init = 0)

# The Huber estimate stays near the uncontaminated mean, unlike mean(y).
coef(m)
#>   theta_1 
#> 0.1567699 
mean(y)
#> [1] 1.082826
```

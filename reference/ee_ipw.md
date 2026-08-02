# Estimating equations for inverse probability weighting (IPW)

Estimates the average causal effect using IPW with a logistic propensity
score model.

## Usage

``` r
ee_ipw(theta, y, A, W, truncate = NULL, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length `3 + b`, where `b` is the number of
  propensity score model parameters.

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n binary treatment indicators (0/1).

- W:

  Numeric n-by-b design matrix for the propensity score model.

- truncate:

  Optional length-2 numeric vector `c(lower, upper)` to clip estimated
  propensity scores. Bounds must be in ascending order
  (`lower <= upper`). Default `NULL`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A `(3+b)`-by-n matrix of estimating equation contributions, with the
first three rows named `ACE`, `E[Y^1]`, and `E[Y^0]` and the propensity
score rows named `W_1` through `W_b` for the columns of `W`.

## Examples

``` r
# A binary treatment, two confounders, and a continuous outcome whose true
# average causal effect is 1.5.
set.seed(42)
n <- 1000
W1 <- rnorm(n)
W2 <- rbinom(n, 1, 0.4)
A <- rbinom(n, 1, inverse_logit(-0.5 + 0.5 * W1 + 0.3 * W2))
Y <- 2 + 1.5 * A + W1 - 0.5 * W2 + rnorm(n)

W_ps <- cbind(1, W1, W2) # Propensity score design matrix

psi <- function(theta) ee_ipw(theta, y = Y, A = A, W = W_ps)

# theta holds the average causal effect, the mean under treatment, and the
# mean under no treatment, followed by the three propensity score
# coefficients.
m <- m_estimate(stacked_equations = psi, init = rep(0, 6))
coef(m)[1:3]
#>      ACE   E[Y^1]   E[Y^0] 
#> 1.575945 3.313529 1.737584 
```

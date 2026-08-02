# Estimating equations for instrumental variable (IV) estimation

Estimates the causal effect using the usual IV / Wald estimator. The
parameter of interest is the additive effect of treatment A on outcome Y
leveraging instrument Z.

## Usage

``` r
ee_iv_causal(theta, y, A, Z, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length 2: the causal effect and the mean of the
  instrument.

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n observed treatment values.

- Z:

  Numeric vector of n binary instrument values (0/1).

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A 2-by-n matrix of estimating equation contributions, with rows named
`causal_effect` and `mean_Z`.

## Examples

``` r
# An unmeasured confounder U biases the association between A and Y, but the
# instrument Z affects Y only through A. The true causal effect is 3.
set.seed(123)
n <- 500
Z <- rbinom(n, 1, 0.5)
U <- rnorm(n)
A <- rbinom(n, 1, inverse_logit(-1 + 3 * Z + U))
Y <- 3 * A - U + rnorm(n, sd = 0.5)

psi <- function(theta) ee_iv_causal(theta, y = Y, A = A, Z = Z)

# theta holds the causal effect followed by the mean of the instrument.
m <- m_estimate(stacked_equations = psi, init = c(0, 0.5))
coef(m)
#> causal_effect        mean_Z 
#>      2.867646      0.470000 
```

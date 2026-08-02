# Estimating equations for g-estimation of structural nested mean models

Estimates the parameters of a structural nested mean model via
g-estimation. Supports both inefficient (X = NULL) and efficient (X
provided) g-estimators, and both linear and log-linear (Poisson)
structural mean models.

## Usage

``` r
ee_gestimation_snmm(
  theta,
  y,
  A,
  W,
  V,
  X = NULL,
  model = "linear",
  weights = NULL
)
```

## Arguments

- theta:

  Numeric vector. For the inefficient g-estimator, length is `b + c`
  (SMM parameters + PS model parameters). For the efficient g-estimator,
  length is `b + c + d` (SMM + PS + outcome model parameters).

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n binary treatment indicators (0/1).

- W:

  Numeric n-by-c design matrix for the propensity score model.

- V:

  Numeric n-by-b design matrix for the structural mean model. Should NOT
  include A itself.

- X:

  Optional n-by-d design matrix for the outcome model (efficient
  g-estimator). Default `NULL` (inefficient g-estimator).

- model:

  Character string: `"linear"` or `"poisson"`. Default `"linear"`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A matrix of estimating equation contributions. The structural mean model
rows are named `SNM phi_0` through `SNM phi_(b-1)`, matching the
zero-based subscripts the literature gives those parameters. The
propensity score rows are named `W_1` through `W_c`, and, for the
efficient g-estimator, the outcome model rows are named `X_1` through
`X_d`.

## Examples

``` r
# A confounded binary treatment whose true effect on the outcome is -2.
set.seed(42)
n <- 500
W <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, 0.25 + 0.5 * W)
Y <- 5 + 2 * W - 2 * A + rnorm(n)

W_ps <- cbind(1, W) # Propensity score design matrix

# An intercept-only structural mean model gives a single causal contrast.
# Build it with rep() so the column has one entry per observation.
V <- cbind(rep(1, n))

psi <- function(theta) {
  ee_gestimation_snmm(theta, y = Y, A = A, W = W_ps, V = V, model = "linear")
}

# theta holds the structural mean model coefficient, which is the causal
# effect, followed by the two propensity score coefficients.
m <- m_estimate(stacked_equations = psi, init = rep(0, 3))
coef(m)
#> SNM phi_0       W_1       W_2 
#> -2.066667 -1.024504  2.234017 
```

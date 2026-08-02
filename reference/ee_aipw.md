# Estimating equations for augmented inverse probability weighting (AIPW)

Estimates the average causal effect using AIPW, which combines a
propensity score model and an outcome model for doubly-robust
estimation.

## Usage

``` r
ee_aipw(theta, y, A, W, X, X1, X0, truncate = NULL, force_continuous = FALSE)
```

## Arguments

- theta:

  Numeric vector of length `3 + b + c`, where `b` is the number of
  propensity score model parameters and `c` is the number of outcome
  model parameters.

- y:

  Numeric vector of n observed outcomes.

- A:

  Numeric vector of n binary treatment indicators (0/1).

- W:

  Numeric n-by-b design matrix for the propensity score model.

- X:

  Numeric n-by-c design matrix for the outcome model.

- X1:

  Numeric n-by-c design matrix under A=1 for all units.

- X0:

  Numeric n-by-c design matrix under A=0 for all units.

- truncate:

  Optional length-2 numeric vector `c(lower, upper)` to clip propensity
  scores. Bounds must be in ascending order (`lower <= upper`). Default
  `NULL`.

- force_continuous:

  Logical. Force linear regression for outcome model? Default `FALSE`.

## Value

A `(3+b+c)`-by-n matrix of estimating equation contributions, with the
first three rows named `ACE`, `E[Y^1]`, and `E[Y^0]`, the propensity
score rows named `W_1` through `W_b`, and the outcome model rows named
`X_1` through `X_c`.

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
X <- cbind(1, A, W1, W2) # Outcome model, observed design matrix
X1 <- cbind(1, 1, W1, W2) # Outcome model, everyone treated
X0 <- cbind(1, 0, W1, W2) # Outcome model, everyone untreated

psi <- function(theta) {
  ee_aipw(theta, y = Y, A = A, W = W_ps, X = X, X1 = X1, X0 = X0)
}

# theta holds the average causal effect, the mean under treatment, and the
# mean under no treatment, followed by the three propensity score
# coefficients and the four outcome model coefficients.
m <- m_estimate(stacked_equations = psi, init = rep(0, 10))
summary(m, subset = 1:3)
#> ── MEstimator Results ──────────────────────────────────────────────────────────
#> Observations: 1000
#> Parameters: 10
#> 
#>          Estimate    Std.Err    Z-score    95% LCL    95% UCL    P-value    S-value
#> ACE        1.5750     0.0666    23.6451     1.4445     1.7056     <2e-16   408.1897
#> E[Y^1]     3.3069     0.0595    55.6246     3.1904     3.4235     <2e-16        Inf
#> E[Y^0]     1.7319     0.0542    31.9673     1.6257     1.8381     <2e-16   742.4752
```

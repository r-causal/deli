# Estimating equations for the g-formula (g-computation)

Returns a stacked set of estimating equations for the g-formula. When
`X0 = NULL`, estimates a single causal mean under the plan encoded by
`X1`. When `X0` is provided, estimates the average causal effect
(difference between two plans).

## Usage

``` r
ee_gformula(theta, y, X, X1, X0 = NULL, force_continuous = FALSE)
```

## Arguments

- theta:

  Numeric vector. If `X0 = NULL`, length is `1 + p` (causal mean +
  regression coefficients). If `X0` is provided, length is `3 + p` (ACE,
  mean under X1, mean under X0, regression coefficients).

- y:

  Numeric vector of n observed outcomes.

- X:

  Numeric n-by-p design matrix (observed data).

- X1:

  Numeric n-by-p design matrix under action plan 1.

- X0:

  Optional n-by-p design matrix under action plan 0. Default `NULL`.

- force_continuous:

  Logical. Force linear regression even when `y` is binary? Default
  `FALSE`.

## Value

A matrix of estimating equation contributions. When `X0 = NULL` the
first row is named `causal_mean`; when `X0` is provided the first three
rows are named `ACE`, `E[Y^1]`, and `E[Y^0]`, where 1 and 0 index the
two plans. The outcome model rows are named `X_1` through `X_p` for the
columns of `X`.

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

X <- cbind(1, A, W1, W2) # Observed design matrix
X1 <- cbind(1, 1, W1, W2) # Everyone treated
X0 <- cbind(1, 0, W1, W2) # Everyone untreated

psi <- function(theta) ee_gformula(theta, y = Y, X = X, X1 = X1, X0 = X0)

# theta holds the average causal effect, the mean under treatment, and the
# mean under no treatment, followed by the four outcome model coefficients.
m <- m_estimate(stacked_equations = psi, init = rep(0, 7))
coef(m)[1:3]
#>      ACE   E[Y^1]   E[Y^0] 
#> 1.578994 3.309943 1.730949 
```

# Estimating equation for multinomial logistic regression

Supports unranked categorical outcomes. `y` must be an n-by-k indicator
matrix where the first column is the reference category, and k must be
at least three: a two-level outcome is logistic regression, which
[`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md) fits
with `distribution = "binomial"`. Returns a `(b * (k-1))`-by-n matrix.

## Usage

``` r
ee_mlogit(theta, X, y, weights = NULL, offset = NULL)
```

## Arguments

- theta:

  Numeric vector of length `b * (k-1)`.

- X:

  Numeric n-by-b design matrix.

- y:

  Numeric n-by-k indicator matrix (first column = reference), with k at
  least three.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A `(b * (k-1))`-by-n matrix.

## Examples

``` r
set.seed(123)
n <- 50
W <- rbinom(n, 1, 0.5)
probs <- cbind(0.5 - 0.2 * W, 0.3 + 0.1 * W, 0.2 + 0.1 * W)
y_cat <- sapply(seq_len(n), function(i) sample(1:3, 1, prob = probs[i, ]))

# The outcome is an indicator matrix whose first column is the reference
# category.
y <- cbind(
  as.integer(y_cat == 1),
  as.integer(y_cat == 2),
  as.integer(y_cat == 3)
)
X <- cbind(1, W)

psi <- function(theta) ee_mlogit(theta, X = X, y = y)

# Two columns of X and two non-reference categories give four parameters.
m <- m_estimate(stacked_equations = psi, init = rep(0, 4))
coef(m)
#>   theta_1   theta_2   theta_3   theta_4 
#> -1.098612  1.386294 -1.098612  1.704748 
```

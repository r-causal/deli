# Estimating equation for extended Rogan-Gladen correction

Extended version that conditions sensitivity and specificity on
covariates using logistic regression models.

## Usage

``` r
ee_rogan_gladen_extended(theta, y, y_star, r, X, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length `1 + 2*p`: corrected proportion, then `p`
  sensitivity model parameters, then `p` specificity model parameters.

- y:

  Numeric vector of gold-standard measurements (validation sample where
  `r = 0`).

- y_star:

  Numeric vector of mismeasured outcome values.

- r:

  Numeric indicator: 1 for main study, 0 for validation.

- X:

  Numeric n-by-p design matrix for sensitivity/specificity models.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A `(1+2*p)`-by-n matrix, with rows named `corrected_proportion`, then
`sens_1` through `sens_p` for the sensitivity model, then `spec_1`
through `spec_p` for the specificity model.

## Examples

``` r
# A validation design of the kind ee_rogan_gladen() takes, with sensitivity
# and specificity now modeled by logistic regression. The design matrix here
# is intercept only, so both models estimate a single log-odds.
set.seed(1)
n <- 500
y_true <- rbinom(n, 1, 0.3)
y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
r <- c(rep(0, 200), rep(1, 300))

# The gold standard is observed only in the validation sample, so the main
# study positions carry a 0 placeholder, as on ee_rogan_gladen().
y <- ifelse(r == 0, y_true, 0)
X <- cbind(rep(1, n))

psi <- function(theta) {
  ee_rogan_gladen_extended(theta, y = y, y_star = y_star, r = r, X = X)
}

m <- m_estimate(stacked_equations = psi, init = c(0.5, 1, 1))

# Corrected prevalence, then the sensitivity and specificity intercepts
coef(m)
#> corrected_proportion               sens_1               spec_1 
#>            0.2963931            2.0614230            1.7176515 
```

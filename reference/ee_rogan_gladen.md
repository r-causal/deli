# Estimating equation for Rogan-Gladen correction

Corrects for mismeasured binary outcomes using external validation data
to estimate sensitivity and specificity.

## Usage

``` r
ee_rogan_gladen(theta, y, y_star, r, weights = NULL)
```

## Arguments

- theta:

  Numeric vector of length 4: corrected proportion, naive proportion,
  sensitivity, specificity.

- y:

  Numeric vector of gold-standard measurements (only available in
  external validation sample where `r = 0`).

- y_star:

  Numeric vector of mismeasured outcome values (all observations).

- r:

  Numeric indicator: 1 for main study data, 0 for external validation.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

## Value

A 4-by-n matrix, with rows named `corrected_proportion`,
`naive_proportion`, `sensitivity`, and `specificity`.

## Examples

``` r
# A main study measures a binary outcome with an imperfect test. An external
# validation study measures both the test and the gold standard, and so
# informs the sensitivity and specificity used to correct the main study.
set.seed(2)
n_main <- 500
n_validation <- 400
n <- n_main + n_validation
y_true <- rbinom(n, 1, 0.25)
y_star <- ifelse(y_true == 1, rbinom(n, 1, 0.9), 1 - rbinom(n, 1, 0.85))
r <- c(rep(1, n_main), rep(0, n_validation))

# The gold standard is observed only in the validation sample, so the main
# study positions carry a 0 placeholder. It never reaches an estimate: the
# sensitivity and specificity equations are multiplied by (1 - r).
y <- ifelse(r == 0, y_true, 0)

psi <- function(theta) {
  ee_rogan_gladen(theta, y = y, y_star = y_star, r = r)
}

m <- m_estimate(
  stacked_equations = psi,
  init = c(0.5, 0.5, 0.75, 0.75)
)

# Corrected prevalence, naive prevalence, sensitivity, specificity
coef(m)
#> corrected_proportion     naive_proportion          sensitivity 
#>            0.2586982            0.3440000            0.9134615 
#>          specificity 
#>            0.8547297 
```

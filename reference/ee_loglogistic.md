# Estimating equation for 4-parameter log-logistic dose-response model

Implements the 4-parameter log-logistic model: \$\$R_i = \theta_0 +
\frac{\theta_m - \theta_0}{1 + \exp\[\theta_s (\log D_i - \log
\theta\_{50})\]}\$\$

## Usage

``` r
ee_loglogistic(theta, dose, response, loss = NULL, k = NULL)
```

## Arguments

- theta:

  Numeric vector of length 4: lower limit, upper limit, ED50, steepness.

- dose:

  Numeric vector of n dose values.

- response:

  Numeric vector of n response values.

- loss:

  Optional character string for robust loss function.

- k:

  Optional numeric tuning parameter for robust loss.

## Value

A 4-by-n matrix, with rows named `lower`, `upper`, `ed50`, and
`steepness`.

## See also

[`ee_loglogistic_ed()`](https://r-causal.github.io/deli/reference/ee_loglogistic_ed.md)
for the effective dose at a given level, which is stacked with this
equation to give it a sandwich standard error.

## Examples

``` r
# Four-parameter log-logistic dose-response for the ryegrass data.
psi <- function(theta) {
  ee_loglogistic(theta, dose = inderjit$dose, response = inderjit$response)
}

# The default rootSolve solver does not converge here, so nleqslv is used.
m <- m_estimate(
  stacked_equations = psi,
  init = c(0.2, 8, 2, 1),
  solver = "nleqslv"
)

# Lower limit, upper limit, ED50, and steepness
coef(m)
#>     lower     upper      ed50 steepness 
#>  0.481410  7.792962  3.057955  2.982229 
```

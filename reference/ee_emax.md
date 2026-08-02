# Estimating equation for E-max dose-response model

Implements the hyperbolic E-max (Hill) model: \$\$R_i = \theta_0 +
\frac{\theta\_{max} D_i}{\theta\_{50} + D_i}\$\$

## Usage

``` r
ee_emax(theta, dose, response, loss = NULL, k = NULL)
```

## Arguments

- theta:

  Numeric vector of length 3: zero-dose response (`e0`), maximum change
  in response (`emax`), ED50. `emax` is the change the response
  approaches as the dose grows without bound, not the response itself,
  so the asymptote is `e0 + emax`.

- dose:

  Numeric vector of n dose values.

- response:

  Numeric vector of n response values.

- loss:

  Optional character string for robust loss function. Default `NULL` (no
  robust loss). See
  [`robust_loss_functions()`](https://r-causal.github.io/deli/reference/robust_loss_functions.md).

- k:

  Optional numeric tuning parameter for robust loss.

## Value

A 3-by-n matrix, with rows named `e0`, `emax`, and `ed50`.

## See also

[`ee_emax_ed()`](https://r-causal.github.io/deli/reference/ee_emax_ed.md)
for the effective dose at a given level, which is stacked with this
equation to give it a sandwich standard error.

## Examples

``` r
# Dose-response of a herbicide on ryegrass root length. The response falls
# with dose, so the maximum change in response is initialized negative.
psi <- function(theta) {
  ee_emax(theta, dose = inderjit$dose, response = inderjit$response)
}

m <- m_estimate(stacked_equations = psi, init = c(8, -8, 2))

# Zero-dose response, maximum change in response, and ED50
coef(m)
#>        e0      emax      ed50 
#>  8.215129 -9.820041  4.574524 
```

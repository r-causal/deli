# Estimating equation for generalized linear models

Returns a p-by-n matrix of estimating equation contributions for a GLM
with specified distribution and link function: \$\$\psi_i(\theta) =
\\Y_i - g^{-1}(X_i^T \theta)\\ \frac{D(\theta)}{v(\theta)} X_i\$\$

## Usage

``` r
ee_glm(
  theta,
  X,
  y,
  distribution,
  link,
  hyperparameter = NULL,
  weights = NULL,
  offset = NULL
)
```

## Arguments

- theta:

  Numeric vector of length p.

- X:

  Numeric n-by-p design matrix.

- y:

  Numeric vector of n observed outcome values.

- distribution:

  Character string: `"normal"` (or `"gaussian"`), `"binomial"` (or
  `"bernoulli"`, or `"bin"`), `"poisson"`, `"gamma"`,
  `"negative_binomial"` (or `"nb"`), `"inverse_normal"` (or
  `"inverse_gaussian"`), `"tweedie"`. Each alias is one Python
  delicatessen accepts, and this equation accepts the same set.

- link:

  Character string: `"identity"`, `"log"`, `"logit"` (or `"logistic"`),
  `"probit"`, `"cauchit"` (or `"cauchy"`), `"loglog"`, `"cloglog"`,
  `"inverse"`, `"sqrt"` (or `"square_root"`).

- hyperparameter:

  Numeric scalar power \\p\\ for the variance function of the
  `"tweedie"` distribution. Must be non-negative. Ignored by every other
  distribution. Default `NULL`.

- weights:

  Optional numeric vector of n weights. Default `NULL`.

- offset:

  Optional numeric vector of n offsets. Default `NULL`.

## Value

A p-by-n matrix. For the gamma and negative binomial distributions the
matrix has `ncol(X) + 1` rows, the final row being the nuisance equation
for that distribution's extra parameter (the gamma shape or the negative
binomial dispersion). Those two distributions name their rows `X_1`
through `X_p` for the columns of `X`, followed by `log_shape` or
`log_dispersion`. Every other distribution estimates one coefficient per
design column and leaves its rows unnamed.

## Gamma

The gamma distribution estimates one additional shape parameter beyond
the regression coefficients, so `theta` has length `ncol(X) + 1`. The
final element is the log of the shape \\\alpha\\ (the reciprocal of the
GLM dispersion \\\phi = 1 / \alpha\\), and the variance function is
\\\mu^2\\. An extra nuisance estimating equation for \\\log(\alpha)\\ is
appended: \$\$ \left(1 - \frac{Y_i}{\mu_i}\right) +
\log\left(\frac{\alpha Y_i}{\mu_i}\right) - \psi(\alpha) \$\$ where
\\\psi\\ is the digamma function. The returned matrix therefore has
`ncol(X) + 1` rows.

## Negative binomial

The negative binomial distribution estimates one additional dispersion
parameter beyond the regression coefficients, so `theta` has length
`ncol(X) + 1`. The final element is the log of the dispersion
\\\alpha\\, and the variance function is \\\mu + \alpha \mu^2\\. An
extra nuisance estimating equation for \\\log(\alpha)\\ is appended so
that the uncertainty in the dispersion is carried honestly through the
sandwich variance: \$\$ -\alpha^{-2}\psi\left(Y_i + \alpha^{-1}\right) +
\alpha^{-2}\psi\left(\alpha^{-1}\right) + \frac{Y_i}{\alpha^2 \mu_i +
\alpha} - \frac{\frac{\alpha \mu_i}{\alpha \mu_i + 1} +
\log\left(\frac{1}{\alpha \mu_i + 1}\right)}{\alpha^2} \$\$ where
\\\psi\\ is the digamma function. The returned matrix therefore has
`ncol(X) + 1` rows.

## Tweedie

The tweedie distribution uses the variance function \\v(\mu) = \mu^p\\,
where the power \\p\\ is the fixed `hyperparameter` rather than an
estimated parameter. No nuisance estimating equation is appended, so
`theta` has length `ncol(X)` and the returned matrix has `ncol(X)` rows.
The power recovers familiar special cases: \\p = 1\\ gives the Poisson
variance \\\mu\\ and \\p = 2\\ gives the gamma variance \\\mu^2\\, so
those choices share the corresponding beta score equations. The
`hyperparameter` must be non-negative.

## Examples

``` r
# Negative binomial GLM with a log link estimates the regression
# coefficients plus one log-dispersion parameter, so init has ncol(X) + 1
# entries.
set.seed(1)
n <- 200
X <- cbind(1, rnorm(n), rnorm(n))
mu <- exp(0.5 + 0.5 * X[, 2] - 0.3 * X[, 3])
y <- rnbinom(n, size = 2, mu = mu)

psi <- function(theta) {
  ee_glm(theta, X = X, y = y, distribution = "negative_binomial",
         link = "log")
}

m <- m_estimate(stacked_equations = psi, init = c(0, 0, 0, 0))
coef(m)
#>            X_1            X_2            X_3 log_dispersion 
#>      0.3852993      0.6128560     -0.4713772     -0.7962357 

# Tweedie GLM with a log link. The power p is a fixed hyperparameter, not an
# estimated parameter, so init has ncol(X) entries. Here p = 1.5 sits in the
# compound Poisson-gamma regime, appropriate for non-negative data with a
# point mass at zero.
set.seed(2)
mu <- exp(0.5 + 0.5 * X[, 2] - 0.3 * X[, 3])
y_tw <- rpois(n, lambda = mu)

psi_tw <- function(theta) {
  ee_glm(theta, X = X, y = y_tw, distribution = "tweedie", link = "log",
         hyperparameter = 1.5)
}

m_tw <- m_estimate(stacked_equations = psi_tw, init = c(0, 0, 0))
coef(m_tw)
#>    theta_1    theta_2    theta_3 
#>  0.5140316  0.4225417 -0.2271453 
```

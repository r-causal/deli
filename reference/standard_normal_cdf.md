# Standard normal CDF

Evaluates the cumulative distribution function of the standard normal
distribution.

This is the equivalent of
[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html) at the default
mean and standard deviation, and returns identical values for numeric
input, except that it carries derivatives. Exact differentiation
(`deriv_method = "exact"`) propagates a tangent alongside each value,
and `standard_normal_cdf()` recognizes a tangent-carrying argument and
applies the analytic rule itself, the standard normal density.
[`pnorm()`](https://rdrr.io/r/stats/Normal.html) hands its argument to
compiled code without dispatching, and errors on such an argument. Use
`standard_normal_cdf()` inside estimating equations and inside
transforms passed to
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md),
and [`pnorm()`](https://rdrr.io/r/stats/Normal.html) for ordinary
numeric work.

## Usage

``` r
standard_normal_cdf(x)
```

## Arguments

- x:

  Numeric value or vector of quantiles.

## Value

Numeric CDF values.

## Exact differentiation

`deriv_method = "exact"` is forward-mode automatic differentiation: it
replaces each value with an object carrying both the value and its
derivative. deli supports those objects through S3 methods: the `Ops`,
`Math`, and `Summary` group generics, plus non-group methods such as
`[`, `%*%`, [`t()`](https://rdrr.io/r/base/t.html),
[`c()`](https://rdrr.io/r/base/c.html), and
[`mean()`](https://rdrr.io/r/base/mean.html). `standard_normal_cdf()`,
[`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md),
[`deli_polygamma()`](https://r-causal.github.io/deli/reference/deli_polygamma.md),
and
[`deli_digamma()`](https://r-causal.github.io/deli/reference/deli_digamma.md)
recognize a tangent-carrying argument themselves and apply their own
analytic rule. Support within the group generics is partial; see
[`vignette("getting-started")`](https://r-causal.github.io/deli/articles/getting-started.md)
for the operations deli differentiates.

[`plogis()`](https://rdrr.io/r/stats/Logistic.html),
[`qlogis()`](https://rdrr.io/r/stats/Logistic.html),
[`pnorm()`](https://rdrr.io/r/stats/Normal.html),
[`dnorm()`](https://rdrr.io/r/stats/Normal.html), and
[`psigamma()`](https://rdrr.io/r/base/Special.html) take none of these
paths. Each hands its argument straight to compiled code through
[`.Call()`](https://rdrr.io/r/base/CallExternal.html) or
[`.Internal()`](https://rdrr.io/r/base/Internal.html) without
dispatching, so the tangent reaches C code that requires a plain number.
deli catches the resulting failure and raises its own error, naming the
function that stopped the computation and the deli function to write in
its place. The same applies to every other distribution function in
`stats`, such as [`qnorm()`](https://rdrr.io/r/stats/Normal.html), which
is named in the error even though deli exports no counterpart for it.

Each deli utility below returns the same values as its base R
counterpart for numeric input. What separates them is whether the
counterpart survives exact mode:

|  |  |  |
|----|----|----|
| deli function | base R counterpart | base R under `deriv_method = "exact"` |
| [`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md) | [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html) | errors |
| [`logit()`](https://r-causal.github.io/deli/reference/logit.md) | [`stats::qlogis()`](https://rdrr.io/r/stats/Logistic.html) | errors |
| `standard_normal_cdf()` | [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html) | errors |
| [`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md) | [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) | errors |
| [`deli_polygamma()`](https://r-causal.github.io/deli/reference/deli_polygamma.md) | [`base::psigamma()`](https://rdrr.io/r/base/Special.html) | errors |
| [`deli_digamma()`](https://r-causal.github.io/deli/reference/deli_digamma.md) | [`base::digamma()`](https://rdrr.io/r/base/Special.html) | works |
| [`identity_transform()`](https://r-causal.github.io/deli/reference/identity_transform.md) | [`base::identity()`](https://rdrr.io/r/base/identity.html) | works |

For the first five rows, use the deli function inside estimating
equations and inside transforms passed to
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md),
and the base R function everywhere else: simulating data, post-fit
display, plain numeric work. For the last two rows the base R function
is usable everywhere, because
[`digamma()`](https://rdrr.io/r/base/Special.html) is a `Math` group
member with a tangent rule and
[`identity()`](https://rdrr.io/r/base/identity.html) passes its argument
through untouched.

The polygamma row is the one place where the arguments do not line up.
`deli_polygamma(n, x)` takes the derivative order first and
`psigamma(x, deriv = n)` takes it second, so a positional substitution
between the two computes a different quantity and raises no error.

## See also

[`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html), the base R
equivalent for ordinary numeric work, and
[`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md)
for the density.

## Examples

``` r
standard_normal_cdf(0)
#> [1] 0.5
standard_normal_cdf(c(-1.96, 0, 1.96))
#> [1] 0.0249979 0.5000000 0.9750021

# The same values as the base R counterpart
all.equal(standard_normal_cdf(c(-1.96, 0, 1.96)), pnorm(c(-1.96, 0, 1.96)))
#> [1] TRUE

# A probit regression, whose mean model is the standard normal CDF
m <- m_estimate(
  vs ~ mpg,
  data = mtcars,
  .ee = ee_glm,
  distribution = "binomial",
  link = "probit"
)

# Variance of the fitted probability at mpg = 20. Writing
# `pnorm(theta[1] + theta[2] * 20)` here instead would error, because exact
# differentiation hands the transform a tangent-carrying argument.
delta_method(
  m,
  transform = function(theta) standard_normal_cdf(theta[1] + theta[2] * 20),
  deriv_method = "exact"
)
#>            [,1]
#> [1,] 0.01286948
```

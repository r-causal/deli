# Polygamma function

Evaluates the `n`th derivative of the digamma function.

This is the equivalent of
[`base::psigamma()`](https://rdrr.io/r/base/Special.html) and returns
identical values for numeric input, except that it carries derivatives.
Exact differentiation (`deriv_method = "exact"`) propagates a tangent
alongside each value, and `deli_polygamma()` recognizes a
tangent-carrying argument and applies the analytic rule itself, the
polygamma function of order `n + 1`.
[`psigamma()`](https://rdrr.io/r/base/Special.html) hands its argument
to compiled code without dispatching, and errors on such an argument.
Use `deli_polygamma()` inside estimating equations and inside transforms
passed to
[`delta_method()`](https://r-causal.github.io/deli/reference/delta_method.md),
and [`psigamma()`](https://rdrr.io/r/base/Special.html) for ordinary
numeric work.

**The two take their arguments in opposite orders.**
`deli_polygamma(n, x)` takes the derivative order first, matching
Python's `scipy.special.polygamma(n, x)`; `psigamma(x, deriv = n)` takes
it second. `deli_polygamma(1, 3)` is `psigamma(3, deriv = 1)`, not
`psigamma(1, deriv = 3)`. A positional substitution between the two
computes a different quantity and raises no error.

## Usage

``` r
deli_polygamma(n, x)
```

## Arguments

- n:

  Integer order of the derivative of the digamma function.

- x:

  Numeric value or vector.

## Value

Numeric polygamma values.

## Exact differentiation

`deriv_method = "exact"` is forward-mode automatic differentiation: it
replaces each value with an object carrying both the value and its
derivative. deli supports those objects through S3 methods: the `Ops`,
`Math`, and `Summary` group generics, plus non-group methods such as
`[`, `%*%`, [`t()`](https://rdrr.io/r/base/t.html),
[`c()`](https://rdrr.io/r/base/c.html), and
[`mean()`](https://rdrr.io/r/base/mean.html).
[`standard_normal_cdf()`](https://r-causal.github.io/deli/reference/standard_normal_cdf.md),
[`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md),
`deli_polygamma()`, and
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
| [`standard_normal_cdf()`](https://r-causal.github.io/deli/reference/standard_normal_cdf.md) | [`stats::pnorm()`](https://rdrr.io/r/stats/Normal.html) | errors |
| [`standard_normal_pdf()`](https://r-causal.github.io/deli/reference/standard_normal_pdf.md) | [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) | errors |
| `deli_polygamma()` | [`base::psigamma()`](https://rdrr.io/r/base/Special.html) | errors |
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

[`base::psigamma()`](https://rdrr.io/r/base/Special.html), the base R
equivalent for ordinary numeric work, which takes its arguments in the
other order, and
[`deli_digamma()`](https://r-causal.github.io/deli/reference/deli_digamma.md)
for the order-zero case.

## Examples

``` r
deli_polygamma(0, 1)
#> [1] -0.5772157
deli_polygamma(1, c(1, 2, 5))
#> [1] 1.6449341 0.6449341 0.2213230

# The same values as the base R counterpart, with the arguments swapped
all.equal(deli_polygamma(1, c(1, 2, 5)), psigamma(c(1, 2, 5), deriv = 1))
#> [1] TRUE

# A gamma GLM, whose last parameter is the log of the shape
m <- m_estimate(
  mpg ~ wt,
  data = mtcars,
  .ee = ee_glm,
  distribution = "gamma",
  link = "log",
  init = c(3, 0, 0)
)

# Variance of the trigamma function at the estimated shape. Writing
# `psigamma(exp(theta[3]), deriv = 1)` here instead would error, because
# exact differentiation hands the transform a tangent-carrying argument.
delta_method(
  m,
  transform = function(theta) deli_polygamma(1, exp(theta[3])),
  deriv_method = "exact"
)
#>             [,1]
#> [1,] 1.61749e-05
```

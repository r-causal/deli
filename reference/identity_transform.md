# Identity transformation

Returns the input unchanged. Used as a no-op transformation in contexts
that accept an arbitrary transformation function.

This is [`base::identity()`](https://rdrr.io/r/base/identity.html) under
a name that does not mask it. The two are interchangeable everywhere,
including under exact differentiation (`deriv_method = "exact"`), since
returning the argument untouched also returns any tangent it carries.
`identity_transform()` exists so that code translated from Python
delicatessen, where the function is called
[`identity()`](https://rdrr.io/r/base/identity.html), can keep its
shape. In R code, prefer
[`identity()`](https://rdrr.io/r/base/identity.html).

## Usage

``` r
identity_transform(value)
```

## Arguments

- value:

  Any value.

## Value

The input `value`, unchanged.

## See also

[`base::identity()`](https://rdrr.io/r/base/identity.html), which does
the same thing and is the idiomatic choice in R.
[`inverse_logit()`](https://r-causal.github.io/deli/reference/inverse_logit.md)
has the full list of deli utilities and their base R counterparts.

## Examples

``` r
identity_transform(42)
#> [1] 42
identity_transform(c(1, 2, 3))
#> [1] 1 2 3

# The same value as base R, and the same object
identical(identity_transform(c(1, 2, 3)), identity(c(1, 2, 3)))
#> [1] TRUE
```

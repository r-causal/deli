# Forward-mode automatic differentiation

Computes the exact Jacobian of a function using forward-mode automatic
differentiation via primal-tangent pairs. Uses one forward pass per
parameter direction with scalar tangents, which correctly handles
interactions between parameters and data vectors.

## Usage

``` r
auto_differentiation(theta, f)
```

## Arguments

- theta:

  Numeric vector of parameter values at which to evaluate the Jacobian.

- f:

  A function that takes a numeric vector and returns a numeric vector
  (or scalar). The function may use standard arithmetic operators and
  math functions (`+`, `-`, `*`, `/`, `^`, `exp`, `log`, `sqrt`, `sin`,
  `cos`, etc.). The matrix product `%*%`, the transpose
  [`t()`](https://rdrr.io/r/base/t.html), the concatenation
  [`c()`](https://rdrr.io/r/base/c.html),
  [`mean()`](https://rdrr.io/r/base/mean.html), and two-dimensional
  indexing differentiate exactly in any function, because each is a
  registered S3 method. The reshaping, binding, and reduction operations
  [`matrix()`](https://rdrr.io/r/base/matrix.html),
  [`rbind()`](https://rdrr.io/r/base/cbind.html),
  [`cbind()`](https://rdrr.io/r/base/cbind.html),
  [`rep()`](https://rdrr.io/r/base/rep.html), and
  [`colSums()`](https://rdrr.io/r/base/colSums.html) differentiate
  exactly for functions evaluated within the package (such as the
  built-in estimating equations), because those masked forms are only in
  scope there; a user-defined function that calls them from the global
  environment reaches the base versions, and the differentiation aborts
  rather than returning a silent approximation. Coercion cannot preserve
  a derivative: [`as.numeric()`](https://rdrr.io/r/base/numeric.html),
  [`as.double()`](https://rdrr.io/r/base/double.html),
  [`as.integer()`](https://rdrr.io/r/base/integer.html),
  [`as.character()`](https://rdrr.io/r/base/character.html),
  [`as.complex()`](https://rdrr.io/r/base/complex.html),
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html), and
  [`as.vector()`](https://rdrr.io/r/base/vector.html) asked by `mode`
  for one of those types or for `"raw"` each return plain values, which
  have nowhere to carry one, so each aborts on a tangent-carrying value
  rather than dropping the tangents silently. The `"pairlist"` and
  `"expression"` modes abort for the neighboring reason: each hands back
  the primal and the tangent as two ordinary elements of a container,
  which no later operation reads as a derivative.
  [`as.vector()`](https://rdrr.io/r/base/vector.html) with its default
  `mode = "any"` keeps the tangents and drops the dimensions of the
  value it is given, as it does for a plain matrix, so it flattens a
  matrix-shaped result such as `X %*% theta` to a vector.
  [`c()`](https://rdrr.io/r/base/c.html) flattens one too.
  [`drop()`](https://rdrr.io/r/base/drop.html) is not dispatched at all,
  so it is masked within the package as the reshaping operations above
  are: inside the package it removes the extents of length one from a
  tangent-carrying value, and a function that calls it from the global
  environment reaches base R's, which hands the value back whole.
  Flatten with [`as.vector()`](https://rdrr.io/r/base/vector.html) or
  [`c()`](https://rdrr.io/r/base/c.html) there.
  [`as.logical()`](https://rdrr.io/r/base/logical.html) is the
  exception, because a logical coercion is a step function whose
  derivative is zero almost everywhere, so it returns the logical its
  payload coerces to; `as.vector(x, "logical")` asks for the same
  coercion and answers the same way. `log(x, base)` differentiates with
  respect to `x` only; the `base` argument is treated as a constant, and
  a `base` that itself carries a tangent (a value derived from `theta`)
  is not supported and aborts rather than dropping the base's
  contribution. [`median()`](https://rdrr.io/r/stats/median.html),
  [`quantile()`](https://rdrr.io/r/stats/quantile.html), and
  `mean(x, trim)` select among the order statistics of their argument,
  whose derivative is that of whichever order statistic the current
  values select rather than that of the population quantity, so each
  aborts as well.

## Value

A matrix where element `[i, j]` is the partial derivative of the `i`-th
output with respect to the `j`-th parameter.

## Conditions

The aborts raised under exact differentiation carry condition classes,
which are the first in the package: `deli_exact_tangent_lost` when
derivative information is gone or a coercion would discard it,
`deli_exact_unsupported_function` when a function reached with a
tangent-carrying argument has no rule under exact differentiation, and
`deli_exact_unsupported_shape` when a result keeps its tangents but
arrives in a container that summing the estimating equations has no
reduction for.

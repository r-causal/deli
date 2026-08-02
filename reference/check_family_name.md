# Check that a distribution or link is a single string

The generalized linear estimating equations read `distribution` and
`link` by name, and the helpers that read them dispatch through `if` on
a comparison with a single value. A `NULL` or a vector longer than one
never reaches the unsupported-name diagnostic those helpers raise: it
fails the `if` first, as base R's `argument is of length zero` or
`the condition has length > 1`, reported against a branch the caller
never wrote.
[`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md)
partitions `theta` on a comparison of its own before either helper is
reached, so the name is judged where the caller supplied it instead.

## Usage

``` r
check_family_name(value, arg, call = rlang::caller_env())
```

## Arguments

- value:

  The distribution or link supplied by the caller.

- arg:

  The argument name, used in the error message.

- call:

  The frame to report the refusal against, which is the estimating
  equation the caller wrote.

## Value

Invisible `NULL`. Raises an error if the value is not a single
non-missing string.

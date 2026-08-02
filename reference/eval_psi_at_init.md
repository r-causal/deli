# Evaluate and validate the estimating function at the initial values

Performs the single evaluation of the estimating function at `init` that
[`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
needs, and validates the return with
[`check_psi_at_init()`](https://r-causal.github.io/deli/reference/check_psi_at_init.md).

## Usage

``` r
eval_psi_at_init(
  psi,
  init,
  allow_over_identification = FALSE,
  error_call = NULL
)
```

## Arguments

- psi:

  The estimating-function closure.

- init:

  The initial parameter vector.

- allow_over_identification:

  Logical. When `TRUE` (the GMM case), the estimating function may
  return more equations than parameters, so only a shortfall is
  rejected. Default `FALSE`.

- error_call:

  The frame to report a failure at the caller's own `init` against,
  which is the
  [`estimate()`](https://r-causal.github.io/deli/reference/estimate.md)
  method the caller reached. At an `init` the formula interface
  generated, the entry point recorded on the closure is preferred to it.

## Value

The value of `psi(init)`. Raises an error if it is not a valid
estimating-function return. A failure reframed as a problem with the
automatic length carries the class `deli_formula_auto_init_error`, so a
caller can recognize it without matching the message.

## Details

When the formula interface generated `init` itself it records the vector
on the closure as the `deli_auto_init` attribute, and this reframes a
failure at exactly those starting values against the automatic `init`,
which the user never chose. Some estimating equations append parameters
beyond the regression coefficients (for example
[`ee_glm()`](https://r-causal.github.io/deli/reference/ee_glm.md) with
`"gamma"` or `"negative_binomial"`, which add a scale or dispersion
parameter), so the automatic `init` is one short and the estimating
function either fails outright or returns the wrong number of rows.
Where the equation is one the formula interface recognizes, the reframed
message names the parameter the automatic length leaves out; otherwise
how much the message can claim about the length depends on the failure.
A wrong-shaped return is a length mismatch by definition and is
described as one, an error the estimating function raised is described
as one where the error itself reads as a length problem, and an error
that reads as anything else is left described as itself, so that a user
whose equation failed for a reason of its own is not sent looking for a
length that was never wrong.

Only those two failures are reframed. A `NULL` or non-numeric return
keeps its own message, because an `init` one element short cannot cause
either. A non-finite return keeps its own message when the number of
rows fits the automatic length, and is reframed when it does not,
because an estimating function reading a parameter the automatic `init`
does not reach returns `NA`s and the wrong number of rows together, and
the row count is the more accurate of the two. Every failure raised
after this point, in the solver or the sandwich components, keeps its
own message as well.

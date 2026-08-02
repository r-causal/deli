# Check the estimating-function return at the initial values

Evaluates the value of `stacked_equations(init)` and rejects returns
that would otherwise fail deep inside the solver with an unhelpful
message: a `NULL` or non-numeric return, a number of estimating
equations that cannot be solved against the number of parameters, or a
non-finite value at the starting values. This mirrors the up-front
validation Python Delicatessen performs before solving.

## Usage

``` r
check_psi_at_init(
  vals,
  init,
  allow_over_identification = FALSE,
  error_call = NULL
)
```

## Arguments

- vals:

  The value of `stacked_equations(init)`, evaluated once by the caller.

- init:

  The initial parameter vector.

- allow_over_identification:

  Logical. When `TRUE`, the number of estimating equations may exceed
  the number of parameters (the GMM case) and only a shortfall is
  rejected. Default `FALSE`.

- error_call:

  The frame to report the error against. `NULL` reports no call, which
  is what a caller with no frame worth naming leaves.

## Value

Invisible `NULL`. Raises an error if the return is invalid. A mismatch
between the number of estimating equations and the number of parameters
carries the class `deli_psi_shape_error`, which is the one failure here
that an automatically generated `init` can explain. The GMM shortfall
carries the number of moment conditions as the `n_moments` field of that
condition as well.

## Details

The order the checks are judged in is load-bearing, because a single
return can fail more than one of them and only the first is reported.
The shape is two branches with different guards, an exact match under
M-estimation and a shortfall under GMM, and each holds its position for
its own reason. Both sit behind the `NULL` and numeric checks and ahead
of the finite one, so a return whose shape cannot be solved is reported
in preference to a non-finite value at the starting values, under GMM as
much as under M-estimation. See the comments in the body for why each
check sits where it does.

Every abort here judges the estimating function the caller supplied,
from several frames below the method the caller reached. This frame
names the parameters of an internal helper rather than anything in the
call the user made, so each abort reports the frame its caller passes
instead.

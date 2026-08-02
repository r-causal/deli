# Check that two design matrices have identical dimensions

Validates that a counterfactual or plan design matrix has the same
dimensions as the observed design matrix, mirroring the explicit shape
checks Python Delicatessen performs in the causal estimating equations
before any arithmetic can silently recycle.

## Usage

``` r
check_design_dims_match(x, y, x_arg, y_arg)
```

## Arguments

- x:

  The reference design matrix.

- y:

  The design matrix to compare against `x`.

- x_arg:

  The reference argument name, used in the error message.

- y_arg:

  The compared argument name, used in the error message.

## Value

Invisible `NULL`. Raises an error if the dimensions differ.

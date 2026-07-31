#' Compute the bread matrix
#'
#' Computes the bread matrix for the empirical sandwich variance estimator.
#' The bread is the negative Jacobian of the summed estimating equations. It is
#' returned unscaled: the callers that assemble a sandwich divide it by the
#' number of observations, and the meat with it.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions.
#' @param theta Numeric vector of parameter estimates.
#' @param deriv_method Character string for the derivative method. One of
#'   `"capprox"` (central), `"fapprox"` (forward), or `"bapprox"` (backward).
#' @param dx Numeric step size (default `1e-9`). The step is absolute, floored
#'   at the floating-point resolution of each estimate; see
#'   [approx_differentiation()].
#' @param summed_equations A function of `theta` returning the length-p vector of
#'   row sums of `stacked_equations` at `theta`, or `NULL` (default) to derive
#'   that reduction from the full p-by-n return. See [compute_sandwich()] for
#'   what a supplied reduction saves and what it must satisfy.
#'
#'   Anything that is neither `NULL` nor a function, and any return at `theta`
#'   that is not numeric or holds fewer values than there are parameters, raises
#'   an error carrying the class `deli_summed_equations_error`. The values
#'   themselves are taken on trust: this function evaluates the estimating
#'   equations nowhere, so it has nothing to compare them against, and a
#'   reduction that sums some other system returns the Jacobian of that other
#'   system with nothing to say so.
#'
#' @returns The negated Jacobian of the summed estimating equations, with one
#'   row per estimating equation and one column per parameter. That is p-by-p
#'   for an M-estimation system, which has one equation per parameter, and
#'   n_eqs-by-p for an over-identified GMM system, whose rectangular bread
#'   [build_sandwich()] pseudo-inverts. No scaling is applied here; the division
#'   by n that puts the bread on the mean scale belongs to the callers that
#'   assemble a sandwich, [compute_sandwich()] and [estimate()].
#'
#'   A bread holding `NA` is returned as it stands, alongside a warning carrying
#'   the class `deli_bread_na`. What to do about it is the caller's, and the two
#'   callers differ: a fit records no variance and says so, while
#'   [compute_sandwich()] has nothing but the matrix to return and fails.
#'
#' @keywords internal
compute_bread <- function(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-9,
  summed_equations = NULL
) {
  deriv_method <- check_deriv_method(deriv_method)
  check_summed_function(summed_equations)

  # Sum the estimating equations across observations
  derived_summed_ee <- function(input_theta) {
    ef <- stacked_equations(input_theta)
    # Handle PrimalTangent returns (from autodiff-compatible EE functions)
    if (is_pt(ef)) {
      return(sum(ef))
    }
    # A tangent-carrying matrix (the usual p-by-n return of a built-in EE under
    # the exact pass): sum across observations within each equation, carrying
    # the tangents in parallel. Checked before the list branch because a
    # PrimalTangentArray is structurally a list whose `[[` yields scalar pairs.
    #
    # A dimensionless tangent array (a psi built with `c()`, or any vector
    # return) sums every element into one value, whereas the list branch below
    # sums within each element. The two agree only when there is a single
    # equation, which is the only case either shape can reach: `estimate()`
    # treats a vector return as one 1-by-n equation, so a p > 1 fit whose psi
    # returns a vector fails on the resulting 1-by-p Jacobian regardless.
    if (is_pt_array(ef)) {
      if (is.null(dim(ef$primal))) {
        return(primal_tangent(sum(ef$primal), sum(ef$tangent)))
      }
      # `rowSums()` keeps the row labels a psi assigned to name its parameters,
      # and the bread reports none: `estimate()` reads parameter names off the
      # plain numeric evaluation it makes at the solved values, never off a
      # differentiated one. Stripping them here rather than leaving it to the
      # single `as.vector()` that reads the Jacobian column keeps the invariant
      # off any one chokepoint.
      return(primal_tangent_array(
        unname(rowSums(ef$primal)),
        unname(rowSums(ef$tangent))
      ))
    }
    if (is.list(ef) && any(vapply(ef, is_pt, logical(1)))) {
      # A plain list holding one value per equation, from an `lapply()` or
      # `sapply()` inside the psi. `c()` returns a PrimalTangentArray and so takes
      # the branch above instead.
      #
      # Every element is reduced with `sum()`, which each tangent surface has a
      # method for, so an equation that arrives as a scalar pair, as a tangent
      # array, or as a plain numeric constant all sum across observations here.
      # A scalar pair anywhere in the list is what says the list is per-equation:
      # reading the first element alone decided that for the whole list, so a psi
      # whose theta-free equation or whose `c()`-built equation came first was
      # refused the reduction this branch carries out.
      #
      # The names an `lapply()` over named equations records ride the reduction
      # into the Jacobian column and out into the bread's dimnames, so they are
      # dropped for the same reason the row sums above are.
      check_equation_list(ef, length(theta))
      return(unname(lapply(ef, sum)))
    }
    # Under the exact pass, any other list shape is one this step cannot sum.
    # `base::rbind()` on a single tangent-carrying value, which is what a psi
    # written in the global environment reaches, builds a 1-by-2 list matrix
    # holding one slot per cell. Neither cell holds a pair, which is why it misses
    # the per-equation branch above; its `dim` is why it would go on to reach
    # `rowSums()` below, which fails with the opaque `'x' must be numeric or
    # complex`. Route it to a diagnostic instead.
    #
    # The two ways a list can arrive here are different failures. A list whose
    # elements still carry tangents, such as an `lapply()` returning one tangent
    # array per equation, has lost nothing and only its container shape is
    # unsupported, so it gets its own condition rather than being told the
    # derivatives are gone.
    if (deriv_method == "exact" && is.list(ef)) {
      if (any(vapply(ef, is_tangent_container, logical(1)))) {
        pt_unsupported_shape_abort()
      }
      pt_tangent_lost_abort()
    }
    if (is.null(dim(ef))) {
      # A plain list under a finite-difference method is the per-equation return
      # the exact pass reduces above, evaluated at plain numbers, so it takes the
      # same reduction. Nothing succeeds under the exact pass that would fail
      # without it, and handing this shape to `sum()` reported base R's
      # `invalid 'type' (list) of argument` from a frame no deli rule runs in.
      if (is.list(ef)) {
        check_equation_list(ef, length(theta))
        return(unname(vapply(ef, sum, numeric(1))))
      }
      return(sum(ef))
    }
    unname(rowSums(ef))
  }

  # The reduction above is what the derivative is taken of, and a caller who
  # supplied one of their own replaces the whole of it. Nothing else changes:
  # the derived reduction is the same closure it always was, so a call that
  # supplies nothing differentiates exactly what it differentiated before.
  summed_ee <- summed_equations %||% derived_summed_ee

  # One plain evaluation of a supplied reduction, read before one or two per
  # parameter are spent differentiating it. This function evaluates nothing else,
  # so the shape is the whole of what it can judge; whether the reduction sums
  # the estimating functions it was passed with is [compute_sandwich()]'s
  # reading, and a caller that comes here directly is trusted on it.
  if (!is.null(summed_equations)) {
    check_summed_return(
      summed_equations(theta),
      length(theta),
      exact = FALSE
    )
  }

  if (deriv_method == "exact") {
    bread_matrix <- auto_differentiation(theta, summed_ee)
  } else {
    bread_matrix <- approx_differentiation(
      func = summed_ee,
      theta = theta,
      method = deriv_method,
      dx = dx
    )
  }

  if (anyNA(bread_matrix)) {
    cli::cli_warn(
      c(
        "!" = "The bread matrix contains NA values, so it cannot be inverted.",
        "i" = "The variance will not be calculated."
      ),
      class = "deli_bread_na"
    )
  }

  -1 * bread_matrix
}

# Judge a per-equation list before anything is built from it. Two rules, both of
# them about what one element stands for, and both raised here so that the bread
# and the meat take exactly the same lists.
#
# The count. Each element is reduced with a single `sum()` for the bread, so each
# contributes exactly one row, and an element holding a block of several
# equations still reduces to one value. Such an element silently stands in for
# every equation it holds, and the bread comes back with fewer rows than the
# system has, which is a shape nothing downstream can tell from a correct one.
# A list return holds one element per equation, and an M-estimation system has
# one equation per parameter, so counting elements against parameters catches it.
#
# The lengths. The equations of one system are evaluated at one sample, so every
# element holds one value per observation and all of them are the same length.
# An element of another length is an equation built from the wrong observations,
# and neither half of the sandwich reports it on its own: the reduction sums each
# element whatever its length, and `rbind()` recycles a short element up to the
# width of the longest, silently where that width is a multiple of it.
#
# The count is judged first because it is the coarser reading of the same
# mistake, and its message names the block an element stands in for. A block
# whose count happens to fit, two elements for two parameters where one of them
# holds two equations, is left to the lengths: it is twice the width of its
# neighbor.
#
# Neither rule can tell a correctly shaped element from a one-row block by the
# payload alone, since an equation over n observations and a one-row block over
# the same n look alike and both are right.
#' @noRd
check_equation_list <- function(ef, n_params) {
  if (length(ef) != n_params) {
    cli::cli_abort(
      c(
        "Summing the estimating equations received a list of {length(ef)}
         element{?s} for {n_params} parameter{?s}.",
        "i" = "Each element of a list return is one estimating equation, summed
               across observations on its own, so the list holds one element per
               parameter.",
        "i" = "An element holding a block of several equations, such as a matrix
               of rows, sums to a single value and yields a bread with too few
               rows. Return the whole system as one p-by-n value instead, built
               with {.fn rbind} or {.code t()} and arithmetic."
      ),
      class = "deli_exact_unsupported_shape",
      call = NULL
    )
  }
  # `lengths()` reads the `length()` method of each element, which answers for
  # the payload a tangent-carrying value holds rather than for its two slots, so
  # this reads observations under the exact pass as it does under the others.
  widths <- unique(lengths(ef))
  if (length(widths) > 1L) {
    cli::cli_abort(
      c(
        "The estimating equations arrived in a list holding {length(widths)}
         different lengths.",
        "i" = "Each element of a list return is one estimating equation
               evaluated at every observation, so each holds one value per
               observation and all of them are the same length.",
        "i" = "An element of another length is either an equation built from a
               different sample or a block of several equations. Return the
               whole system as one p-by-n value instead, built with {.fn rbind}
               or {.code t()} and arithmetic."
      ),
      class = "deli_exact_unsupported_shape",
      call = NULL
    )
  }
  invisible(ef)
}

# Bind a per-equation list into the p-by-n matrix the meat is built from.
# `compute_bread()` reduces the same list with one `sum()` per element, which is
# the derivative it needs; the meat needs those equations unreduced, so the
# elements are bound into rows here instead. What the two share is the judgment
# above rather than the reduction, which is why it is the whole of this
# function besides the binding.
#' @noRd
stack_equation_list <- function(ef, theta) {
  check_equation_list(ef, length(theta))
  do.call(rbind, ef)
}

#' Compute the meat matrix
#'
#' Computes the meat matrix as the cross-product of the estimating equation
#' evaluations: \eqn{EE \times EE^T}.
#'
#' @param evaluations A p-by-n matrix of estimating equation evaluations, where
#'   p is the number of parameters and n is the number of observations.
#'
#' @returns A p-by-p meat matrix.
#'
#' @keywords internal
compute_meat <- function(evaluations) {
  tcrossprod(evaluations)
}

#' Build the sandwich variance estimator
#'
#' Combines bread and meat matrices into the sandwich:
#' \eqn{B^{-1} M (B^{-1})^T}.
#'
#' @param bread A bread matrix with one row per estimating equation and one
#'   column per parameter. It is p-by-p for an M-estimation system, which has
#'   one equation per parameter, and n_eqs-by-p for an over-identified GMM
#'   system, whose rectangular bread has no inverse and is pseudo-inverted
#'   instead.
#' @param meat An n_eqs-by-n_eqs meat matrix, square in the estimating equations
#'   whether or not the bread is.
#' @param allow_pinv Logical. If `TRUE` (default), uses the pseudo-inverse
#'   when the bread matrix cannot be inverted. When `FALSE`, a bread that has no
#'   inverse raises an error carrying the class `deli_bread_not_invertible`.
#' @param call The frame to report that error against.
#'
#' @returns A p-by-p sandwich covariance matrix, or `NULL` if the bread
#'   contains `NA` values.
#'
#' @keywords internal
build_sandwich <- function(
  bread,
  meat,
  allow_pinv = TRUE,
  call = rlang::caller_env()
) {
  # An NA bread is returned as a fit with no variance rather than as a failure,
  # which is the state check_estimated() names and the accessors report. The
  # entry point that has nothing but this matrix to hand back refuses the return
  # before the assembly is reached; see compute_sandwich().
  if (anyNA(bread)) {
    return(NULL)
  }

  # Invert the bread matrix
  if (allow_pinv) {
    # Use pseudo-inverse for robustness
    bread_inv <- tryCatch(
      solve(bread),
      error = function(e) {
        rlang::check_installed(
          "MASS",
          reason = "for pseudo-inverse when the bread matrix is singular."
        )
        MASS::ginv(bread)
      }
    )
  } else {
    # The shape is read ahead of the solve because it is not a question the
    # solve can answer usefully: base R reports a rectangular bread as
    # `'a' (2 x 1) must be square`, naming an argument of its own, and an
    # over-identified system has a remedy of its own to be told about. Whether a
    # square bread has an inverse is left to the solve, which is the whole of
    # what that question means here.
    check_bread_square(bread, call = call)
    bread_inv <- tryCatch(
      solve(bread),
      error = function(e) abort_bread_not_invertible(bread, call = call)
    )
  }

  # Sandwich: B^{-1} M (B^{-1})^T
  bread_inv %*% meat %*% t(bread_inv)
}

# ---- the bread that cannot be inverted ---------------------------------------
# Two conditions report a bread the sandwich cannot be assembled from, and both
# carry a class so that a caller or a test can match on the class rather than on
# the prose, as the solver warnings in `R/estimate.R` do.
#
#   deli_bread_na
#     The bread holds `NA`, so no inverse of it exists. A warning, raised by
#     compute_bread(), because what follows from it depends on the caller: a fit
#     records no variance and carries on, and compute_sandwich() converts it into
#     the error below. See compute_bread().
#
#   deli_bread_not_invertible
#     The bread cannot be inverted and the caller refused the pseudo-inverse, or
#     asked for a matrix and there is none to give. Raised by
#     check_bread_square() for a rectangular bread under `allow_pinv = FALSE`,
#     by abort_bread_not_invertible() where the solve of a square one fails
#     under the same setting, and by compute_sandwich() for a bread holding
#     `NA`.
#
#   deli_meat_not_invertible
#     The covariance of the moment conditions cannot be inverted and the caller
#     refused the pseudo-inverse, so the two-step GMM weight update has no
#     weight matrix to take. Raised by abort_meat_not_invertible(), from the
#     update in estimate_gmm_estimator().

#' Refuse a bread that has no inverse at any rank
#'
#' A rectangular bread is the over-identified GMM system, whose bread has no
#' inverse however well conditioned it is, and base R reported it as
#' `'a' (2 x 1) must be square`, naming an argument of its own. It is refused
#' here with the reason it has no inverse and the setting that would accept one.
#'
#' This is the whole of what is read ahead of the solve. Whether a square bread
#' can be inverted is what the solve answers, and answering it another way is
#' answering a different question: the rank `qr()` returns is read at a
#' tolerance of `1e-7` relative to the largest pivot, so reading it here refused
#' every bread whose condition number runs past roughly `1e7`, including the
#' ones base R inverts and the caller has a more discriminating gate of their
#' own for. deli's contract is solve-or-fail, and `solve()` fails at a
#' reciprocal condition number near the resolution of a double. See
#' `abort_bread_not_invertible()`, which is what the failure reaches.
#'
#' @param bread The bread matrix, with one row per estimating equation and one
#'   column per parameter.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_bread_not_invertible` where the bread is not square.
#' @noRd
check_bread_square <- function(bread, call = rlang::caller_env()) {
  n_eqs <- nrow(bread)
  n_params <- ncol(bread)
  if (n_eqs == n_params) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "The bread matrix has {n_eqs} row{?s} and {n_params} column{?s},
             so it has no inverse.",
      "i" = "An over-identified system has more estimating equations than
             parameters, and its rectangular bread is pseudo-inverted rather
             than solved.",
      "i" = "Set {.code allow_pinv = TRUE} to build the variance from the
             Moore-Penrose pseudo-inverse, which is what an over-identified
             fit reports."
    ),
    class = "deli_bread_not_invertible",
    call = call
  )
}

#' Refuse a bread the solve could not invert
#'
#' The counterpart of `abort_meat_not_invertible()` for the other half of the
#' sandwich, and raised the same way: where the solve fails rather than ahead of
#' it. That is what makes the refusal cover both ways a bread defeats the
#' inverse, and only those ways. Linearly dependent estimating equations are the
#' usual one and the factorization names the directions. The other is a bread
#' whose columns are independent to the factorization's tolerance and whose
#' scales differ by more than the reciprocal condition number `solve()` accepts,
#' which a rank reading cannot see at all.
#'
#' What a rank reading ahead of the solve sees instead is every bread whose
#' condition number runs past roughly `1e7`, since that is where the `1e-7`
#' tolerance `qr()` applies falls. A finite-difference Jacobian of a system with
#' a dependent direction lands there routinely: the round-off in the difference
#' quotient perturbs the dependent row, so the matrix reads as rank deficient and
#' inverts without complaint. Refusing it took the reading out of the caller's
#' hands, and the caller's own is the more discriminating one. A dependence
#' confined to a nuisance block leaves the effects a fit reports identified to
#' every digit they had, and a fit is still told about the dependence: see
#' `not_identified()`, which reads the same matrix at the same tolerance and
#' warns rather than refusing to return a variance.
#'
#' The rank is read here, on the path where the matrix has already been found to
#' have no inverse, because it is what says which directions were lost and it
#' costs one factorization of a matrix the call is about to abandon.
#'
#' @param bread The bread matrix the solve could not invert.
#' @param call The frame to report the error against.
#'
#' @returns Nothing; it throws.
#' @noRd
abort_bread_not_invertible <- function(bread, call = rlang::caller_env()) {
  n_params <- ncol(bread)
  if (!all(is.finite(bread))) {
    detail <- c(
      "i" = "It holds values that are not finite, so it has no inverse and no
             factorization to read one from."
    )
  } else {
    factored <- qr(bread)
    if (factored$rank < n_params) {
      # The directions are named by the pivoting `qr()` does, which drops one
      # column per lost direction and reports which. Those are the parameters
      # the remaining ones already account for rather than the ones at fault:
      # a dependence is a property of a set, and which member of the set is
      # named is the factorization's choice. The positions are written out as
      # strings, and the count stated with `qty()`, for the reason
      # `warn_dependent_moments()` writes its own out: cli would otherwise read
      # a single number as the quantity to pluralize on.
      dependent <- as.character(sort(factored$pivot[-seq_len(factored$rank)]))
      detail <- c(
        "i" = "Its rank is {factored$rank} of {n_params}, so it is singular: at
               least one direction in the parameter space leaves the mean
               estimating equations unchanged, and the parameters along it
               cannot be told apart.",
        "i" = "{cli::qty(dependent)}Parameter{?s} {dependent} {?is/are}
               accounted for by the others."
      )
    } else {
      detail <- c(
        "i" = "It factors at full rank, so no direction in the parameter space
               is accounted for by the others, and what defeats the inverse is
               the range of scale across them."
      )
    }
  }
  cli::cli_abort(
    c(
      "!" = "The bread matrix cannot be inverted, and {.arg allow_pinv} is
             {.code FALSE}.",
      detail,
      "i" = "Set {.code allow_pinv = TRUE} to build the variance from the
             Moore-Penrose pseudo-inverse, or look for a redundant parameter, a
             design whose columns are linearly dependent, and an estimating
             equation whose scale stands far from the rest."
    ),
    class = "deli_bread_not_invertible",
    call = call
  )
}

#' Refuse a moment covariance that has no inverse
#'
#' The counterpart of `check_bread_invertible()` for the two-step GMM weight
#' update, which sets the weight matrix to the inverse of the moment covariance.
#' `allow_pinv = FALSE` says a matrix with no inverse is to be refused rather
#' than pseudo-inverted, and the update took its inverse with a bare
#' [base::solve()], so a covariance with none surfaced as LAPACK's
#' `system is exactly singular` or `system is computationally singular`, naming
#' neither the setting that produced the refusal nor the conditions at fault.
#'
#' The refusal is raised where the solve fails rather than ahead of it, which is
#' what makes it cover both ways a covariance defeats the inverse. Linearly
#' dependent conditions are the usual one and the factorization names them. The
#' other is a covariance whose columns are independent to the factorization's
#' tolerance and whose scales differ by more than the reciprocal condition
#' number `solve()` accepts; reading the rank ahead of the solve would let that
#' one through to fail as base R.
#'
#' @param meat The moment covariance the update could not invert.
#' @param call The frame to report the error against.
#'
#' @returns Nothing; it throws.
#' @noRd
abort_meat_not_invertible <- function(meat, call = rlang::caller_env()) {
  detail <- if (!all(is.finite(meat))) {
    "It holds values that are not finite, so it has no inverse and no
     factorization to read one from."
  } else {
    factored <- qr(meat)
    n_moments <- ncol(meat)
    if (factored$rank < n_moments) {
      "Its rank is {factored$rank} of {n_moments}, so at least one moment
       condition is accounted for by the others and their covariance has no
       inverse."
    } else {
      "It factors at full rank, so no one condition is accounted for by the
       others, and what defeats the inverse is the range of scale across them."
    }
  }
  cli::cli_abort(
    c(
      "!" = "The moment covariance cannot be inverted, and {.arg allow_pinv} is
             {.code FALSE}.",
      "i" = detail,
      "i" = "Set {.code allow_pinv = TRUE} to build the two-step weight matrix
             from the Moore-Penrose pseudo-inverse, which is what an
             over-identified fit reports, or drop the moment conditions the
             others already account for."
    ),
    class = "deli_meat_not_invertible",
    call = call
  )
}

#' Apply finite-sample correction to the meat matrix
#'
#' Applies the HC1 correction: \eqn{meat \times n / (n - p)}.
#'
#' @param meat A p-by-p meat matrix.
#' @param n Integer number of observations.
#' @param p Integer number of parameters.
#' @param adjustment Character string or `NULL`. Currently only `"HC1"` is
#'   supported.
#'
#' @returns The corrected meat matrix.
#'
#' @keywords internal
finite_sample_correction <- function(meat, n, p, adjustment = NULL) {
  if (is.null(adjustment)) {
    return(meat)
  }

  adj_upper <- toupper(adjustment)
  if (adj_upper == "HC1") {
    if (n <= p) {
      cli::cli_abort(
        "The number of observations ({n}) is not greater than the number of
         parameters ({p}), so the HC1 correction cannot be applied."
      )
    }
    meat * n / (n - p)
  } else {
    cli::cli_abort(
      c(
        "The finite-sample correction {.val {adjustment}} is not available.",
        "i" = "Supported options: {.val NULL}, {.val HC1}."
      )
    )
  }
}

# ---- the reduction the bread is differentiated from --------------------------
# Four readings of `summed_equations` and the switch that governs the last of
# them, in the order a call reaches them, and one family class across all of
# them. [MEstimator()] and [GMMEstimator()] carry the pair as properties and
# reach the same readings, the first two at construction and the third where the
# fit has the evaluation to make it.
#
#   check_summed_function()
#     The argument itself. Nothing evaluates it until the bread does, so an
#     argument that is not a function surfaced from wherever it was first called
#     as base R's `could not find function`, naming an internal formal.
#
#   check_summed_switch()
#     The switch. It is a property rather than an argument on the fit path, so
#     it is read once at construction rather than at every use, and a value that
#     is not a single `TRUE` or `FALSE` would otherwise be collapsed to one by
#     `isTRUE()` with nothing said.
#
#   check_summed_return()
#     The shape of one evaluation. A return that is not numeric failed inside the
#     difference quotient as `non-numeric argument to binary operator`, and a
#     return of the wrong length built a Jacobian with too few rows, which has no
#     inverse at any rank and nothing to say why.
#
#   check_summed_agreement()
#     Whether the values are the ones the meat is built from, which only
#     compute_sandwich() has the evaluation to judge.

#' Check that a supplied reduction is a function
#'
#' @param summed_equations The argument as the caller passed it.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_summed_equations_error` for anything that is neither `NULL` nor a
#'   function.
#' @noRd
check_summed_function <- function(
  summed_equations,
  call = rlang::caller_env()
) {
  if (is.null(summed_equations) || is.function(summed_equations)) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "{.arg summed_equations} must be a function or {.val NULL}, not
             {.obj_type_friendly {summed_equations}}.",
      "i" = "It takes the parameter vector and returns the sum of each
             estimating equation across the observations, which is what the
             bread is differentiated from."
    ),
    class = "deli_summed_equations_error",
    call = call
  )
}

#' Check that the agreement switch is a single logical value
#'
#' @param check_summed_equations The argument as the caller passed it.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_summed_equations_error` for anything that is not a single `TRUE` or
#'   `FALSE`.
#' @noRd
check_summed_switch <- function(
  check_summed_equations,
  call = rlang::caller_env()
) {
  usable <- is.logical(check_summed_equations) &&
    length(check_summed_equations) == 1L &&
    !is.na(check_summed_equations)
  if (usable) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "{.arg check_summed_equations} must be {.code TRUE} or
             {.code FALSE}.",
      "i" = "It says whether {.arg summed_equations} is compared against the
             estimating functions it is paired with, which is the reading that
             tells a reduction of this system from a reduction of another one."
    ),
    class = "deli_summed_equations_error",
    call = call
  )
}

#' Check the shape of one evaluation of a supplied reduction
#'
#' The two callers know different amounts about the count. [compute_sandwich()]
#' has the estimating functions in hand and so knows exactly how many values the
#' reduction holds, one per row of their return. [compute_bread()] evaluates
#' nothing else and knows only that a system has at least one estimating equation
#' per parameter, since a Jacobian with fewer rows than columns has no inverse at
#' any rank; more of them is the over-identified system, whose rectangular bread
#' is pseudo-inverted.
#'
#' @param summed The value of the reduction at `theta`.
#' @param n_values The number of values the return is judged against.
#' @param exact Whether that number is the count it must hold, or the fewest it
#'   may hold.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the class
#'   `deli_summed_equations_error` for a return the bread cannot be built from.
#' @noRd
check_summed_return <- function(
  summed,
  n_values,
  exact = TRUE,
  call = rlang::caller_env()
) {
  if (!is.numeric(summed)) {
    cli::cli_abort(
      c(
        "!" = "{.arg summed_equations} must return a numeric vector at
               {.arg theta}, not {.obj_type_friendly {summed}}.",
        "i" = "The bread is the Jacobian of what it returns, so the return holds
               one number per estimating equation."
      ),
      class = "deli_summed_equations_error",
      call = call
    )
  }
  n_summed <- length(summed)
  if (if (exact) n_summed == n_values else n_summed >= n_values) {
    return(invisible(NULL))
  }
  detail <- if (exact) {
    "{.arg stacked_equations} returns {n_values} estimating equation{?s} at
     {.arg theta}, and the reduction holds the sum of each of them across the
     observations."
  } else {
    "A system has at least one estimating equation for each of its {n_values}
     parameter{?s}, and the reduction holds the sum of each of them across the
     observations. Fewer leaves the bread with more columns than rows, so it has
     no inverse at any rank."
  }
  cli::cli_abort(
    c(
      "!" = "{.arg summed_equations} returned {n_summed} value{?s} at
             {.arg theta}.",
      "i" = detail
    ),
    class = "deli_summed_equations_error",
    call = call
  )
}

#' Check that a supplied reduction sums the estimating functions it is paired
#' with
#'
#' The correctness hazard `summed_equations` introduces is a reduction that does
#' not sum the estimating functions the meat is built from. The bread would then
#' be the Jacobian of one system and the meat the cross-product of another, and
#' the returned matrix would carry the shape of a covariance without the
#' meaning. Nothing downstream can tell such a matrix from a correct one.
#'
#' The check is close to free, which is why it is on by default. The full
#' evaluation is already in hand for the meat, so one `rowSums()` over it is
#' O(p n), which the call is already paying, rather than the O(p^2 n) the
#' argument exists to avoid. It costs one further call to the reduction.
#'
#' The comparison is relative, with an absolute floor of one. Two summations of
#' the same n values in different orders differ by rounding, which is relative to
#' the size of what was summed, so a fixed absolute tolerance would refuse a
#' correct reduction on a large sample. The floor is there because the row sums
#' vanish at a root, which is the point this function is usually called at, and a
#' purely relative comparison has no scale to stand on there.
#'
#' That the comparison is made at one point, and that the point is a root, is
#' also the limit of what it can see. Both quantities are at rounding there, so a
#' reduction that is a multiple of the right one agrees with it and is accepted,
#' and the bread comes back as that multiple of the right bread. What the check
#' catches is a reduction of some other system, which does not vanish where this
#' one does. Catching a rescaling would take an evaluation of the estimating
#' functions somewhere other than the root, which is the whole cost the argument
#' exists to avoid.
#'
#' @param summed The value of the supplied reduction at `theta`.
#' @param from_matrix The row sums of the full evaluation at `theta`.
#' @param call The frame to report the error against.
#'
#' @returns Invisible `NULL`. Raises an error carrying the classes
#'   `deli_summed_equations_disagree` and `deli_summed_equations_error` where the
#'   two disagree, and the family class alone where the shape is what is wrong.
#' @noRd
check_summed_agreement <- function(
  summed,
  from_matrix,
  call = rlang::caller_env()
) {
  n_eqs <- length(from_matrix)
  check_summed_return(summed, n_eqs, exact = TRUE, call = call)
  summed <- unname(as.numeric(summed))

  # A value that is not a number is read before the comparison rather than
  # through it. The comparison takes the largest disagreement, and `which.max()`
  # passes over what it cannot order: a return that is missing throughout named
  # no equation at all and reported base R's `argument is of length zero` from
  # the test itself, and a return missing in one equation was passed over
  # altogether, differentiated, and surfaced two steps later as a bread that
  # could not be inverted.
  #
  # The values it is compared against are the row sums of a return already judged
  # finite, so nothing here is a reading of the estimating functions: a reduction
  # that is not a number where they are cannot be summing them.
  if (anyNA(summed)) {
    missing <- as.character(which(is.na(summed)))
    cli::cli_abort(
      c(
        "!" = "{.arg summed_equations} returned a value that is not a number at
               {.arg theta}.",
        "i" = "{cli::qty(missing)}Estimating equation{?s} {missing} {?sums/sum}
               to {.val {NA_real_}}, where {.arg stacked_equations} sums to a
               number.",
        "i" = "The bread is differentiated from {.arg summed_equations}, so a
               bread built from this one holds {.val {NA_real_}} throughout and
               no sandwich can be assembled from it."
      ),
      class = c(
        "deli_summed_equations_disagree",
        "deli_summed_equations_error"
      ),
      call = call
    )
  }

  # A comparison of no values has nothing to disagree about, and an ordering of
  # none has no equation to name.
  if (n_eqs == 0L) {
    return(invisible(NULL))
  }
  scale <- pmax(abs(from_matrix), 1)
  relative <- abs(summed - from_matrix) / scale
  # Subtracting one infinity from another is the one remaining way the
  # difference of two numbers is not one, and it takes a row sum that overflowed
  # on both sides. Ordering passes over such a value, which would again leave the
  # largest disagreement naming no equation, so it is read as the disagreement it
  # is.
  relative[is.na(relative)] <- Inf
  worst <- which.max(relative)
  if (relative[[worst]] < 1e-6) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "{.arg summed_equations} does not sum {.arg stacked_equations} at
             {.arg theta}.",
      "i" = "Estimating equation {worst} sums to
             {.val {from_matrix[[worst]]}} across the observations, and
             {.arg summed_equations} reports {.val {summed[[worst]]}}.",
      "i" = "The bread is differentiated from {.arg summed_equations} and the
             meat is built from {.arg stacked_equations}, so the two must be the
             same system."
    ),
    class = c("deli_summed_equations_disagree", "deli_summed_equations_error"),
    call = call
  )
}

#' Compute the empirical sandwich variance estimator
#'
#' Computes the empirical sandwich variance estimator directly from a set of
#' estimating equations and a vector of parameter estimates. Unlike
#' [MEstimator()], this function does not solve for the parameters; it assumes
#' that `theta` is already the root of the estimating equations and only
#' assembles the sandwich covariance at that point.
#'
#' The sandwich is built from a bread matrix and a meat matrix. The bread is the
#' negative Jacobian of the summed estimating equations, and the meat is the
#' cross-product of the equation evaluations. Each is scaled by `1/n` internally,
#' so the returned matrix is on the asymptotic scale (see `Value`). The bread
#' Jacobian is obtained either by finite differences or by forward-mode
#' automatic differentiation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of parameters and n is the number of observations. A list of one
#'   element per equation, each holding that equation's contributions across the
#'   observations, is accepted as well. [estimate()] has no support for that
#'   shape, so an estimating function written in it reaches a variance through
#'   this entry point and a Jacobian through [compute_bread()]. The list form
#'   holds exactly one element per parameter, each of one length, so an
#'   over-identified system reaches this function as a matrix rather than as a
#'   list.
#' @param theta Numeric vector of parameter estimates. This function assumes
#'   `theta` is the root of `stacked_equations`; it does not solve for it.
#' @param deriv_method Character string selecting the method used to build the
#'   bread Jacobian. One of `"capprox"` (central difference), `"fapprox"`
#'   (forward difference), `"bapprox"` (backward difference), or `"exact"`
#'   (forward-mode automatic differentiation). Default `"capprox"`.
#'
#'   This default differs from Python delicatessen, whose `compute_sandwich`
#'   defaults to `"approx"`, a forward difference computed through SciPy's
#'   `approx_fprime`. deli does not replicate the SciPy `"approx"` path; its
#'   `"fapprox"` is the hand-implemented forward difference. Code ported across
#'   the two libraries should set `deriv_method` explicitly rather than rely on
#'   the default.
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. A small value is recommended, since large steps
#'   can give poor approximations. Default `1e-9`. Must be a single positive
#'   finite number, which is checked whichever `deriv_method` is in force. The
#'   step is absolute and is floored at the floating-point resolution of each
#'   estimate, so a large parameter magnitude cannot silently reduce it to
#'   nothing; see [approx_differentiation()].
#' @param allow_pinv Logical. When `TRUE` (default), the Moore-Penrose
#'   pseudo-inverse is used where the bread matrix cannot be solved, which is
#'   the case for the rectangular bread of an over-identified system as well as
#'   for a square one whose solve fails; when `FALSE`, a bread with no inverse
#'   raises an error carrying the class `deli_bread_not_invertible`.
#'
#'   What the two settings choose between is what to do when the solve fails,
#'   and nothing else: a bread the solve inverts is inverted under either. An
#'   ill-conditioned bread that base R returns an inverse for is one of those,
#'   however far its condition number runs, so `allow_pinv = FALSE` is not a
#'   conditioning test and a caller who wants one applies it to the returned
#'   matrix.
#' @param finite_correction Character string or `NULL`. Finite-sample correction
#'   applied to the meat matrix. `NULL` (default) applies no correction; `"HC1"`
#'   rescales the meat by \eqn{n / (n - p)}, where p is the number of parameters.
#' @param summed_equations A function that takes a numeric vector `theta` and
#'   returns the length-p vector of row sums of `stacked_equations` at `theta`,
#'   or `NULL` (default) to derive those sums from the full p-by-n return.
#'
#'   The bread is the Jacobian of the summed estimating equations, so each of the
#'   one or two perturbed evaluations it makes per parameter is reduced to one
#'   value per equation as soon as it is built. Deriving the reduction builds the
#'   whole p-by-n matrix 2p times for arithmetic that is linear in it; an
#'   estimating function whose sums have a closed form, such as the
#'   \eqn{X^T r} of a regression score, can supply them and hand the bread only
#'   what it uses. The meat is unaffected either way, since it needs the
#'   per-observation contributions and takes the one full evaluation it always
#'   took.
#'
#'   Under `deriv_method = "exact"` the reduction is called with a
#'   tangent-carrying `theta`, so it must be written in operations that carry
#'   derivatives: `t(X) %*% r` does, and [base::crossprod()] does not. See
#'   [auto_differentiation()] for which operations carry a tangent and where.
#'
#'   An argument that is neither `NULL` nor a function raises an error carrying
#'   the class `deli_summed_equations_error`, whatever `check_summed_equations`
#'   says. So does a return at `theta` that is not numeric, or one holding fewer
#'   values than there are parameters, both of which [compute_bread()] reads
#'   before it differentiates anything. That the return holds exactly one value
#'   per estimating equation is read by the comparison below, since only a call
#'   that has evaluated the estimating functions knows how many of them there
#'   are.
#' @param check_summed_equations Logical. When `TRUE` (default) and
#'   `summed_equations` was supplied, its value at `theta` is compared against
#'   the row sums of the one full evaluation the meat is built from, and a
#'   disagreement raises an error carrying the class
#'   `deli_summed_equations_disagree`. The comparison costs one call to
#'   `summed_equations` and one reduction of a matrix the call already holds.
#'   Set it to `FALSE` to skip the comparison, which leaves a reduction that sums
#'   some other system to return a matrix with the shape of a covariance and no
#'   claim to be one.
#'
#'   The comparison is made at `theta`, which the caller states is the root, so
#'   both quantities are at rounding there. That is enough to catch a reduction
#'   of some other system and not enough to catch a reduction that is a multiple
#'   of the right one, which agrees at a root and yields that multiple of the
#'   right bread.
#'
#' @returns A p-by-p covariance matrix on the asymptotic scale. The bread and
#'   meat are each divided by n internally, so the returned matrix is the
#'   variance that corresponds to the standard deviation. Dividing it by the
#'   number of observations gives the standard-error-scale variance, whose
#'   square-rooted diagonal is the vector of standard errors.
#'
#'   A covariance matrix is the only thing this function returns. Where the
#'   bread has no inverse, and so no sandwich can be assembled from it, the call
#'   raises an error carrying the class `deli_bread_not_invertible` rather than
#'   returning something that has to be tested for. That covers a bread holding
#'   `NA`, and, under `allow_pinv = FALSE`, a rectangular bread and one the
#'   solve could not invert. [estimate()] makes the other choice from the same
#'   matrices: a fit whose
#'   bread holds `NA` warns and comes back with no variance, since it still
#'   carries the estimates.
#'
#' @seealso [MEstimator()] and [GMMEstimator()], which solve for `theta` and
#'   report this variance internally, and [delta_method()] for the variance of a
#'   transformation of the parameters.
#'
#' @examples
#' # A generic data set for estimating a mean and variance
#' y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
#'
#' # The mean and variance are the roots of ee_mean_variance, so they can be
#' # computed directly rather than solved for
#' theta <- c(mean(y), stats::var(y) * (length(y) - 1) / length(y))
#'
#' # Wrap the built-in estimating equation as a function of theta alone
#' psi <- function(theta) ee_mean_variance(theta, y = y)
#'
#' # compute_sandwich() returns the asymptotic-scale variance, so dividing by
#' # n puts it on the standard-error scale
#' sandwich <- compute_sandwich(psi, theta = theta) / length(y)
#' sandwich
#'
#' # The diagonal square roots are the standard errors
#' sqrt(diag(sandwich))
#'
#' # The bread only ever needs the sums of the estimating equations, so an
#' # equation whose sums have a closed form can supply them directly
#' summed <- function(theta) {
#'   c(
#'     sum(y) - length(y) * theta[1],
#'     sum((y - theta[1])^2) - length(y) * theta[2]
#'   )
#' }
#' compute_sandwich(psi, theta = theta, summed_equations = summed) / length(y)
#'
#' @references
#' Boos DD, & Stefanski LA. (2013). M-estimation (estimating equations). In
#' Essential Statistical Inference (pp. 297-337). Springer, New York, NY.
#'
#' @export
compute_sandwich <- function(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  finite_correction = NULL,
  summed_equations = NULL,
  check_summed_equations = TRUE
) {
  check_dx(dx)
  # Every failure below judges something the caller passed to this function, and
  # some of them are raised from inside a handler or a helper, where one frame up
  # is no longer this one. The frame is taken once, here, and passed to each of
  # them.
  call <- rlang::current_env()
  deriv_method <- check_deriv_method(deriv_method, call = call)
  check_summed_function(summed_equations, call = call)
  # This evaluates the estimating function once for itself and once or twice per
  # parameter for the bread, so an estimating function that warns raises the same
  # warning several times for one call. See R/conditions.R. The body is short
  # enough to wrap in place, unlike the estimate() methods, which put theirs in a
  # worker.
  without_repeated_warnings({
    # Evaluate estimating equations at theta-hat. The return is judged before
    # anything is built from it, so a return the bread and the meat cannot be
    # assembled from is reported as itself rather than as whatever the assembly
    # goes on to fail at. Over-identification is allowed, because the sandwich
    # this function builds from a rectangular bread is the asymptotic variance a
    # GMM fit reports.
    evald <- stacked_equations(theta)
    # A per-equation list is a documented return of this entry point, the one
    # `compute_bread()` reduces, so it is bound into the matrix the meat is
    # built from before the return is judged. Judging the list itself refused it
    # as a non-numeric return, which meant the bread took a shape the sandwich
    # would not. A classed list is left to the judgment below, where a data
    # frame of contributions is named as a data frame.
    if (is.list(evald) && !is.object(evald)) {
      evald <- stack_equation_list(evald, theta)
    }
    check_psi_at_theta(evald, theta, call = call)
    if (is.null(dim(evald))) {
      n_obs <- length(evald)
      n_params <- 1
      # Reshape to 1-by-n matrix so tcrossprod works correctly in compute_meat
      evald <- matrix(evald, nrow = 1)
    } else {
      n_obs <- ncol(evald)
      n_params <- nrow(evald)
    }

    # A reduction of the estimating functions is judged against the evaluation
    # the meat is built from, before the bread is differentiated from it. See
    # check_summed_agreement() for what the two disagreeing would produce and
    # why reading it here is what makes the check affordable.
    if (!is.null(summed_equations) && isTRUE(check_summed_equations)) {
      check_summed_agreement(
        summed_equations(theta),
        unname(rowSums(evald)),
        call = call
      )
    }

    # Step 1: Bread matrix
    #
    # An NA bread is the one failure the assembly reports rather than raises,
    # because a fit that has the estimates can carry on without a variance. This
    # entry point has nothing else to hand back, so the warning compute_bread()
    # raises is converted into the failure it is here rather than delivered
    # alongside one. Fetching it by class rather than testing the returned matrix
    # is what keeps the report to a single condition: the handler is an exiting
    # one, so the warning never reaches the caller.
    bread <- rlang::try_fetch(
      compute_bread(
        stacked_equations,
        theta,
        deriv_method,
        dx,
        summed_equations = summed_equations
      ),
      deli_bread_na = function(cnd) {
        cli::cli_abort(
          c(
            "!" = "The bread matrix contains NA values, so it has no inverse.",
            "i" = "The sandwich is built from the bread, so there is no
                   covariance matrix to return at {.arg theta}.",
            "i" = "This usually means the estimating functions are not
                   differentiable at {.arg theta}, or that they return a value
                   that is not finite at a perturbed one."
          ),
          class = "deli_bread_not_invertible",
          call = call
        )
      }
    )
    bread <- bread / n_obs

    # Step 2: Meat matrix
    meat <- compute_meat(evald)
    meat <- meat / n_obs
    meat <- finite_sample_correction(meat, n_obs, n_params, finite_correction)

    # Step 3: Build sandwich
    build_sandwich(bread, meat, allow_pinv, call = call)
  })
}

#' Confidence bands for parameter vectors
#'
#' Computes simultaneous confidence bands that provide coverage for the
#' entire parameter vector, adjusting for multiple comparisons. The formula is:
#' \deqn{\hat{\theta} \pm \hat{c}_{\alpha/2} \times \widehat{SE}(\hat{\theta})}
#' where \eqn{\hat{c}} is the adjusted critical value.
#'
#' @param object A fitted `MEstimator` object, or a numeric vector of
#'   parameter estimates.
#' @param alpha Numeric significance level, between 0 and 1. Default `0.05`.
#' @param method Character string. `"supt"` (default) for supremum-t or
#'   `"bonferroni"` for Bonferroni correction.
#' @param n_draws Integer number of MVN draws for the sup-t method.
#'   Default `1e5`. The critical value is a Monte Carlo estimate of a fixed
#'   sup-t quantile, and at `1e5` draws the band half-width varies by roughly
#'   0.2% across independent draws, which is negligible relative to the band
#'   width. Python's estimator method defaults to `1e6`; both defaults target
#'   the same quantity and agree to within Monte Carlo error, so `deli` keeps
#'   the smaller default for faster computation. Set `n_draws = 1e6` to match
#'   Python's estimator default. Note that seeded band values are not
#'   reproducible across the two languages because the MVN samplers differ.
#' @param seed Integer seed for reproducibility. Default `NULL`.
#' @param subset Integer vector of parameter indices to compute bands for.
#'   Default `NULL` (all parameters).
#' @param covariance Numeric covariance matrix (only when `object` is numeric).
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A p-by-2 matrix with columns `"lower"` and `"upper"`. For a fitted
#'   estimator the rows are named for the parameters, as in
#'   [`confint()`][deli-generics]. For a numeric `object` they take their names
#'   from it, when it has any.
#'
#' @examplesIf requireNamespace("MASS", quietly = TRUE)
#' # Two independent samples, each contributing one mean parameter
#' set.seed(42)
#' n <- 200
#' y1 <- rnorm(n, 2)
#' y2 <- rnorm(n, 3)
#'
#' psi <- function(theta) {
#'   rbind(y1 - theta[1], y2 - theta[2])
#' }
#'
#' m <- m_estimate(stacked_equations = psi, init = c(0, 0))
#'
#' # Simultaneous bands cover both means at once, so they are wider than the
#' # pointwise intervals from confint(m)
#' confidence_bands(m, method = "supt", seed = 1)
#'
#' @export
# n_draws defaults to 1e5, not the 1e6 Python uses on its estimator method. The
# sup-t critical value converges to the same quantile either way, and at 1e5 the
# residual Monte Carlo error is negligible relative to the band width while
# running about an order of magnitude faster. Set n_draws = 1e6 for numerical
# agreement with Python under matched settings.
confidence_bands <- new_generic(
  "confidence_bands",
  "object",
  function(
    object,
    alpha = 0.05,
    method = "supt",
    n_draws = 100000L,
    seed = NULL,
    subset = NULL,
    covariance = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

method(confidence_bands, deli_estimator) <- function(
  object,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL,
  subset = NULL,
  covariance = NULL,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)

  if (!is.null(subset)) {
    theta <- object@theta[subset]
    covariance <- object@variance[subset, subset, drop = FALSE]
  } else {
    theta <- object@theta
    covariance <- object@variance
  }

  compute_confidence_bands(
    theta,
    covariance,
    alpha = alpha,
    method = method,
    n_draws = n_draws,
    seed = seed
  )
}

method(confidence_bands, class_numeric) <- function(
  object,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL,
  subset = NULL,
  covariance = NULL,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  if (!is.null(subset)) {
    object <- object[subset]
    if (!is.null(covariance)) {
      covariance <- as.matrix(covariance)[subset, subset, drop = FALSE]
    }
  }

  compute_confidence_bands(
    object,
    covariance,
    alpha = alpha,
    method = method,
    n_draws = n_draws,
    seed = seed
  )
}

#' Compute confidence bands from theta and covariance
#'
#' @param theta Numeric parameter vector.
#' @param covariance Numeric covariance matrix.
#' @param alpha Significance level. Default `0.05`.
#' @param method `"supt"` or `"bonferroni"`. Default `"supt"`.
#' @param n_draws Number of MVN draws for sup-t. Default `1e5`. See
#'   [confidence_bands()] for why this differs from Python's `1e6` estimator
#'   default and how to match it.
#' @param seed RNG seed. Default `NULL`.
#'
#' @returns A p-by-2 matrix with columns `"lower"` and `"upper"`. Rows take
#'   their names from `theta`, when it has any.
#'
#' @examplesIf requireNamespace("MASS", quietly = TRUE)
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' # Bands from the estimates and covariance alone, without the fitted object
#' compute_confidence_bands(coef(fit), covariance = vcov(fit),
#'                          method = "supt", seed = 1)
#'
#' @export
compute_confidence_bands <- function(
  theta,
  covariance,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL
) {
  # as.numeric() strips the names, so they are kept aside for the row labels of
  # the returned bands, which are the same labels confint() puts on its rows.
  param_names <- names(theta)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)
  k <- length(theta)

  # Validate dimensions and alpha before computing, mirroring Python's
  # confidence_bands. Without these, a mismatched theta recycles over the
  # covariance and an out-of-range alpha returns infinite or NaN bands.
  n_rows <- nrow(covariance)
  n_cols <- ncol(covariance)
  if (n_rows != n_cols) {
    cli::cli_abort(
      "{.arg covariance} must be a square matrix, but it has {n_rows} row{?s}
       and {n_cols} column{?s}."
    )
  }
  if (k != n_rows) {
    cli::cli_abort(
      "{.arg theta} and {.arg covariance} must share a dimension, but
       {.arg theta} has length {k} and {.arg covariance} has dimension
       {n_rows}."
    )
  }
  check_alpha(alpha)

  se <- sqrt(diag(covariance))

  method <- tolower(method)
  if (method %in% c("supt", "sup-t")) {
    if (any(se <= 0)) {
      cli::cli_abort(
        "At least one parameter has a standard error of zero or less. The sup-t method cannot be applied."
      )
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }
    rlang::check_installed(
      "MASS",
      reason = "to draw from the multivariate normal for sup-t confidence bands."
    )
    # Draw from multivariate normal
    draws <- MASS::mvrnorm(n = n_draws, mu = rep(0, k), Sigma = covariance)
    if (k == 1) {
      draws <- matrix(draws, ncol = 1)
    }
    # Scale by SE and take absolute value
    scaled <- abs(sweep(draws, 2, se, `/`))
    # Max across parameters for each draw
    ts <- apply(scaled, 1, max)
    # Critical value is (1-alpha) quantile
    cv <- quantile(ts, probs = 1 - alpha, names = FALSE)
  } else if (method == "bonferroni") {
    cv <- qnorm(1 - alpha / (2 * k))
  } else {
    cli::cli_abort(
      "Method {.val {method}} is not supported. Use {.val supt} or {.val bonferroni}."
    )
  }

  cb <- cbind(lower = theta - cv * se, upper = theta + cv * se)
  rownames(cb) <- param_names
  cb
}

#' Delta method for variance of transformed parameters
#'
#' Computes the variance-covariance matrix for a transformation of parameters
#' using the Delta Method:
#' \deqn{Var[g(\theta)] \approx G \Sigma G^T}
#' where \eqn{G} is the Jacobian of \eqn{g} and \eqn{\Sigma} is the
#' covariance matrix of \eqn{\theta}.
#'
#' @param object A fitted `MEstimator` object, or a numeric vector of
#'   parameter estimates.
#' @param transform Function that takes `theta` and returns a numeric vector.
#' @param covariance Numeric covariance matrix (only used when `object` is a
#'   numeric vector).
#' @param deriv_method Character string for the derivative method used to build
#'   the Jacobian of `transform`. One of `"capprox"` (central difference),
#'   `"fapprox"` (forward difference), `"bapprox"` (backward difference), or
#'   `"exact"` (forward-mode automatic differentiation). Default `"capprox"`.
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. Default `1e-9`. Must be a single positive finite
#'   number, which is checked whichever `deriv_method` is in force. The step is
#'   absolute and is floored at the floating-point resolution of each estimate,
#'   so a large parameter magnitude cannot silently reduce it to nothing; see
#'   [approx_differentiation()].
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored. Exact names
#'   matter here because `deriv_method` selects how the Jacobian is built, and a
#'   dropped misspelling would leave the default in place and return a different
#'   variance with nothing to signal the substitution.
#'
#' @returns A covariance matrix for `transform(theta)`.
#'
#' @examples
#' fit <- m_estimate(vs ~ mpg, data = mtcars, .ee = ee_regression,
#'                   model = "logistic")
#'
#' # Variance of the odds ratio for mpg, exponentiating the log-odds coefficient
#' delta_method(fit, transform = function(theta) exp(theta[2]))
#'
#' # The same variance from the estimates and covariance alone
#' delta_method(coef(fit), transform = function(theta) exp(theta[2]),
#'              covariance = vcov(fit))
#'
#' @export
delta_method <- new_generic(
  "delta_method",
  "object",
  function(
    object,
    transform,
    covariance = NULL,
    deriv_method = "capprox",
    dx = 1e-9,
    ...
  ) {
    S7::S7_dispatch()
  }
)

method(delta_method, deli_estimator) <- function(
  object,
  transform,
  covariance = NULL,
  deriv_method = "capprox",
  dx = 1e-9,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)
  # One call differentiates `transform` and so evaluates it several times, once
  # for the shape check and once or twice per parameter for the Jacobian. A
  # transform that warns would otherwise repeat itself for one call. See
  # R/conditions.R. The scope wraps the shared worker rather than this body for
  # the same reason it does in the estimate() methods.
  without_repeated_warnings(
    delta_method_impl(
      object@theta,
      transform,
      object@variance,
      deriv_method,
      dx
    )
  )
}

method(delta_method, class_numeric) <- function(
  object,
  transform,
  covariance = NULL,
  deriv_method = "capprox",
  dx = 1e-9,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  # See the method above for why the worker call sits inside the scope.
  without_repeated_warnings(
    delta_method_impl(object, transform, covariance, deriv_method, dx)
  )
}

#' @noRd
delta_method_impl <- function(
  theta,
  g,
  covariance,
  deriv_method = "capprox",
  dx = 1e-9
) {
  # Every abort below judges an argument the caller passed to one of the two
  # delta_method() methods, and this worker appears in no man page, so each
  # names the frame of the method that was reached rather than a function no
  # user can see. Both methods call this from their own body, so one frame up is
  # the method whichever of them the caller went through.
  call <- rlang::caller_env()

  # Validate covariance presence before any coercion. as.matrix(NULL) otherwise
  # fails with an opaque "'data' must be of a vector type" error before the
  # reachable shape checks below.
  if (is.null(covariance)) {
    cli::cli_abort(
      c(
        "{.arg covariance} is required but was not supplied.",
        "i" = "Pass the covariance matrix of the parameter estimates when
               calling {.fn delta_method} on a numeric vector."
      ),
      call = call
    )
  }
  deriv_method <- check_deriv_method(deriv_method, call = call)
  check_dx(dx, call = call)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  # Validate shapes before computing, mirroring the checks in Python's
  # delta_method. Without these, a bad shape either returns a wrongly shaped
  # matrix or fails with an opaque non-conformable error.
  theta_star <- g(theta)
  star_dim <- dim(theta_star)
  # A plain numeric vector has no dim. R's matrix algebra yields column vectors
  # (n x 1), which are one-dimensional in intent, so those are allowed; a
  # genuine matrix (more than one column) or higher-dimensional array is not.
  if (!is.null(star_dim) && (length(star_dim) > 2 || star_dim[2] > 1)) {
    cli::cli_abort(
      "Output from {.arg transform} must be a one-dimensional vector.",
      call = call
    )
  }
  n_rows <- nrow(covariance)
  n_cols <- ncol(covariance)
  if (n_rows != n_cols) {
    cli::cli_abort(
      "{.arg covariance} must be a square matrix, but it has {n_rows} row{?s}
       and {n_cols} column{?s}.",
      call = call
    )
  }
  v <- length(theta)
  if (n_rows != v) {
    cli::cli_abort(
      "{.arg covariance} and the parameter vector must share a dimension, but
       the parameter vector has length {v} and {.arg covariance} has dimension
       {n_rows}.",
      call = call
    )
  }

  # Delta method: G * Sigma * G^T
  g_prime <- transform_jacobian(theta, g, deriv_method, dx)
  g_prime %*% covariance %*% t(g_prime)
}

#' The Jacobian of a transform at a parameter vector
#'
#' The derivative half of the delta method, separated from the product that
#' follows it because not every caller wants the whole covariance matrix. A
#' caller predicting one quantity per row of a design needs the variance of each
#' row on its own, which is the diagonal of that product, and forming the whole
#' matrix to take its diagonal costs the square of the number of rows in memory.
#' Such a caller takes the Jacobian from here and forms
#' `rowSums((G %*% Sigma) * G)` instead.
#'
#' `deriv_method` is expected already normalized by [check_deriv_method()].
#'
#' @param theta The parameter vector to differentiate at.
#' @param g The transform.
#' @param deriv_method The normalized derivative method.
#' @param dx The finite-difference step.
#' @returns A matrix with one row per element of `g(theta)` and one column per
#'   parameter.
#' @noRd
transform_jacobian <- function(theta, g, deriv_method, dx) {
  # Exact autodiff must receive the raw transform: as.numeric() would strip the
  # tangents and yield a zero Jacobian. extract_tangent_column() already handles
  # every shape the raw transform can return. The finite-difference methods need
  # plain numeric output, so they see the wrapped one.
  if (deriv_method == "exact") {
    return(auto_differentiation(theta, g))
  }
  approx_differentiation(
    func = function(t) as.numeric(g(t)),
    theta = theta,
    method = deriv_method,
    dx = dx
  )
}

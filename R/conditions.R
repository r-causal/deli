#' The condition classes deli raises
#'
#' Where a caller might reasonably want to answer an error or a warning rather
#' than let it stop the work, deli attaches a condition class to it, so that
#' [base::tryCatch()], [base::withCallingHandlers()], or [rlang::try_fetch()]
#' can match on the class rather than on the wording of a message. This page
#' lists those classes and says what each of them reports on.
#'
#' Every class below belongs to an error unless the entry says otherwise. A
#' condition carrying two of them leads with the narrower one, so a handler for
#' either the family or the particular fault catches it.
#'
#' A condition a user's own estimating function raises reaches the caller with
#' whatever class that function gave it. Nothing here renames one.
#'
#' @section What an estimating function returned:
#' - `deli_psi_return_error`: the estimating function returned something no
#'   sandwich can be built from, either at the starting values or at the point
#'   handed to [compute_sandwich()]. Every refusal of a return carries it, so
#'   recognizing one is a match on this rather than on four messages.
#' - `deli_psi_shape_error`: the number of estimating equations returned cannot
#'   be solved against the number of parameters. Leads `deli_psi_return_error`.
#' - `deli_psi_list_unsupported`: the estimating function returned a list of one
#'   element per equation. [compute_bread()] and [compute_sandwich()] take that
#'   shape, and no fit can solve it, so [estimate()] refuses it rather than
#'   failing on the summed equations. Leads `deli_psi_return_error`.
#'
#' @section The formula interface:
#' - `deli_formula_ee_lookup_error`: `.ee` is neither a function nor the name of
#'   one. Either it arrived as some other type entirely, or it arrived as a
#'   character vector naming no one function, because no function of that name
#'   was found or because the vector is not of length one. A handler that means
#'   to read the name it failed to resolve should confirm it has a character
#'   vector first.
#' - `deli_formula_ee_signature_error`: the arguments of `.ee` leave something
#'   the interface fills nowhere to go: `theta`, the model matrix it passes as
#'   `X`, the response it passes by position, or the offset an `offset()` term
#'   in the formula supplies. It also covers an equation that writes formals
#'   after its `...`, which can be filled by name and by nothing else, so the
#'   positional response would fall into the dots and leave them unfilled.
#' - `deli_formula_ee_argument_error`: something in the `...` the caller wrote
#'   cannot be forwarded to `.ee`. On its own it is a name matching no argument
#'   of `.ee` exactly, which is what keeps an abbreviation from reaching an
#'   argument the caller did not mean.
#' - `deli_formula_ee_unnamed_argument`: an argument was forwarded with no name,
#'   and so by position, into whichever argument the response did not take.
#'   Leads `deli_formula_ee_argument_error`.
#' - `deli_formula_ee_reserved_argument`: a name the caller supplied is one the
#'   interface fills itself, either `theta`, `X`, or the argument the response
#'   is passed to. Leads `deli_formula_ee_argument_error`.
#' - `deli_formula_auto_init_error`: the `init` the formula interface generated
#'   itself, one zero per model-matrix column, does not fit the estimating
#'   equation.
#'
#' @section Solving:
#' - `deli_solver_not_converged`, a warning: the solver stopped without solving
#'   the estimating equations, either because it reported a failure of its own
#'   or because the point it returned does not solve them.
#' - `deli_nested_solver_error`: a [rootSolve::multiroot()] solve was asked for
#'   while another one was running, which it cannot be. An error rather than a
#'   warning, because the solve cannot be attempted at all.
#' - `deli_gmm_moments_rejected`, a warning: the Hansen J-statistic of an
#'   over-identified GMM fit stands so far out against its reference
#'   distribution as to all but rule the fit out, which usually means the moment
#'   conditions cannot all hold at one value of the parameters.
#' - `deli_gmm_moments_dependent`, a warning: the moment conditions are linearly
#'   dependent at the estimated values, so the two-step weight matrix is a
#'   pseudo-inverse and what comes back is the fit of the independent conditions
#'   alone. It also carries the rarer covariance that has no inverse for want of
#'   conditioning rather than of rank, which the message tells apart.
#' - `deli_meat_not_invertible`: the covariance of the moment conditions cannot
#'   be inverted, so the two-step GMM weight update has no weight matrix to
#'   take, and `allow_pinv = FALSE` refused the pseudo-inverse. The counterpart
#'   of `deli_bread_not_invertible` for the other half of the sandwich.
#'
#' @section The sandwich variance:
#' - `deli_bread_na`, a warning: the bread matrix holds `NA`, so no inverse of
#'   it exists. A fit records no variance and carries on; [compute_sandwich()]
#'   has nothing but a matrix to return and converts it into the error below.
#' - `deli_bread_not_invertible`: the bread cannot be inverted and
#'   `allow_pinv = FALSE` refused the pseudo-inverse, or it is rectangular and
#'   has no inverse at any rank, or it holds `NA`.
#'
#' @section Differentiation:
#' - `deli_exact_unsupported_function`: under `deriv_method = "exact"`, a
#'   function reached with a tangent-carrying argument has no rule, either
#'   because it hands its argument to compiled code without dispatching or
#'   because deli declines to differentiate it. [auto_differentiation()] names
#'   the replacements to use instead.
#' - `deli_exact_unsupported_shape`: the estimating equations arrived in a
#'   container the summing step has no rule for, either because the tangents
#'   survived in a shape it cannot reduce or because a per-equation list does
#'   not hold one element per parameter, each of one length. Nothing was lost;
#'   the shape is the problem.
#' - `deli_exact_tangent_lost`: derivative information is gone or would be,
#'   either because a tangent-carrying value was asked to become a plain double
#'   or because a result arrived with no tangent and evidence that one had been
#'   dropped on the way.
#' - `deli_finite_difference_lost`, a warning: a finite-difference step changed
#'   the function by less than the floating-point spacing of its values, so an
#'   entry of the Jacobian carries none of the digits of the values it came
#'   from.
#'
#' @section Parameter names:
#' - `deli_param_name_collision`: filling the unnamed entries of a partly named
#'   parameter vector with positional `theta_1`, `theta_2`, ... labels would
#'   repeat a label another parameter in the same vector already carries.
#'
#' @name deli-conditions
NULL

# ---- repeated warnings -------------------------------------------------------
# One fit evaluates its estimating function many times: once at the initial
# values, once at the values the solver returned, and once or twice per parameter
# while the bread matrix is built, on top of whatever the solver itself does. An
# estimating function that warns therefore raises the same warning several times
# over what the user asked for as a single operation. `ee_percentile()`,
# `ee_positive_mean_deviation()` and the penalized regressions all warn that they
# are not differentiable, so all of them repeat.
#
# `delta_method()` differentiates a user transform and repeats for the same
# reason, once for the shape check and once or twice per parameter. It carries
# the scope too. `confidence_bands()` does not, because it takes only estimates
# and a covariance matrix and evaluates no user function at all, and neither does
# `influence_functions()`, which evaluates the estimating function exactly once.
#
# `without_repeated_warnings()` is the scope that collapses those repeats. Two
# things about where it sits are deliberate.
#
# It sits at the boundary of an operation rather than inside any estimating
# equation. An estimating equation cannot tell how many times it is about to be
# called, and only the operation knows where one logical unit of work begins and
# ends. Putting the scope there also covers user-written and third-party
# estimating equations on the same terms as the ones this package provides, and
# needs no cooperation from any of them.
#
# It leaves a direct call to an estimating equation alone, because a direct call
# is one evaluation and has nothing to collapse. That is also the behavior of
# Python delicatessen, which `tests/testthat/test-ee-percentile.R` pins.

#' Deliver each distinct warning once for the duration of an expression
#'
#' Evaluates `expr` with a calling handler that re-signals the first warning
#' carrying a given combination of condition classes and message, and muffles
#' every later one. The record of what has been delivered lives in an
#' environment created by this call, so it is discarded when the call returns and
#' two operations in sequence each report their own warnings.
#'
#' A calling handler, rather than a collect-and-replay pass, is what keeps the
#' rest of R's warning machinery working unchanged. Calling handlers run before
#' the default handling of a warning, so under `options(warn = 2)` the first
#' warning is converted into an error where it was signaled exactly as it would
#' be without this scope, and an enclosing [base::withCallingHandlers()] or
#' [base::tryCatch()] established by the caller still sees it with its own class
#' and call intact.
#'
#' Nested scopes are governed by the outermost one. The inner handler is the more
#' recently established, so it sees a warning first, records it and lets it
#' through; the outer handler then muffles anything it has already delivered. One
#' operation nested inside another therefore reports each distinct warning once
#' in total, which is the promise the outer operation makes on its own.
#'
#' @param expr An expression to evaluate.
#' @returns The value of `expr`.
#' @noRd
without_repeated_warnings <- function(expr) {
  seen <- new.env(parent = emptyenv())
  withCallingHandlers(
    expr,
    warning = function(w) {
      key <- warning_key(w)
      if (is.null(seen[[key]])) {
        assign(key, TRUE, envir = seen)
        return()
      }
      # cnd_muffle() rather than invokeRestart("muffleWarning"), because a
      # condition that inherits from "warning" but was signaled without that
      # restart would turn the restart call into an error. cnd_muffle() returns
      # `FALSE` instead, leaving the warning to propagate, which is the harmless
      # outcome.
      rlang::cnd_muffle(w)
    }
  )
}

#' Key identifying a warning for de-duplication
#'
#' The class vector as well as the message, because the class is part of what a
#' condition is: a caller matching on one class must be able to see it even when
#' a differently classed warning reached the same wording first. The message as
#' well as the class vector, because one class this package raises is shared by
#' warnings that say different things, and the classes cli attaches to an
#' unclassed warning are shared by all of them.
#'
#' Not the call, although base R's own [base::warnings()] keys on it. The cost
#' is real: two warnings that read alike but were raised from different frames
#' collapse to one, and the second is dropped rather than merged. It is accepted
#' because including the call would separate nothing that matters while
#' destroying the de-duplication itself. No warning this package raises carries
#' a call. All of them come from `cli::cli_warn()`, which records a call only
#' when one is passed through to [rlang::warn()], and no site here passes one.
#' `conditionCall()` is therefore empty for every one of them, so the call
#' distinguishes none of them. Among warnings raised by user code with
#' [base::warning()], the reported call is the call to the frame that raised it,
#' and that is a different frame at each site the function is evaluated from: a
#' fit reports `psi(init)` at the initial values and
#' `stacked_equations(full_theta)` after the solve, and [delta_method()] reports
#' `g(theta)` for its shape check and `g(t)` for its finite differences. The
#' formula interface is worse still, because it builds its estimating function
#' with [base::do.call()], which inlines the current parameter values into the
#' constructed call, so no two evaluations ever report the same one.
#'
#' @param w A warning condition.
#' @returns A single string. It is never empty, so it is always a usable name in
#'   the environment holding the seen-set.
#' @noRd
warning_key <- function(w) {
  paste(c("w", class(w), conditionMessage(w)), collapse = "\r")
}

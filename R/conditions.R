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
#' warning is converted into an error where it was signalled exactly as it would
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
      # condition that inherits from "warning" but was signalled without that
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

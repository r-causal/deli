# Helpers for the tests that read the conditions an operation raised rather than
# assert on one of them. Each was written out again in every test file that
# needed it before moving here: three files for the warning collector, five for
# the call reader, and no file for both.

#' Every warning signaled while an expression runs
#'
#' Collects the conditions and muffles each one, so what is counted is what the
#' operation delivered rather than what reached the console, and so no warning
#' escapes into the test summary. Muffling short-circuits R's default warning
#' handling, which is what keeps the count the same under `options(warn = 2)`,
#' the setting the whole suite runs under.
#'
#' Used where one call raises more than one warning, or where the count itself is
#' the assertion. A single expected warning is better asserted with
#' `expect_warning()`.
#'
#' @param expr The expression to evaluate.
#' @returns A list of the warning conditions, in the order they were raised.
collect_warnings <- function(expr) {
  caught <- list()
  withCallingHandlers(
    expr,
    warning = function(w) {
      caught[[length(caught) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  caught
}

#' The call a condition reports, rendered as one string
#'
#' Used by the tests that pin which frame an abort names. The rendering goes
#' through rlang rather than `deparse()`: this package's aborts are raised
#' against S7 method frames whose calls inline the estimator itself, and base
#' `deparse()` fails outright on an inlined S4 object under R 4.3, the oldest R
#' this package supports. rlang writes any value it cannot render as code as a
#' placeholder instead.
#'
#' Two consequences for what callers can match on. S7 inlines a method's
#' arguments into the call as promises, and rlang renders a promise as
#' `<promise>` whether or not it has been forced, so no argument value ever
#' reaches the string: match on argument names here and read the values off the
#' call itself, where indexing forces them. And a long call comes back as
#' several lines, hence the collapse; rlang indents what it wraps, so the joined
#' string carries a run of spaces wherever a line broke, and a `fixed = TRUE`
#' match has to stay inside one argument rather than span the comma between two.
#'
#' @param cnd The condition to read.
#' @returns The reported call as a single string.
reported_call <- function(cnd) {
  paste(rlang::expr_deparse(conditionCall(cnd)), collapse = " ")
}

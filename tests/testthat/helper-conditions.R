# Helper for the tests that read the warnings an operation raised rather than
# assert on one of them. Three test files were carrying the same collector, so it
# lives here instead.

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

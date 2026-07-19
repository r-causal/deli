#' Convert between survival analysis measures
#'
#' Converts a survival probability (and optionally a hazard) to other survival
#' analysis metrics.
#'
#' @param survival Numeric survival probability or vector of probabilities.
#' @param hazard Numeric hazard value or vector. Required for `"hazard"` and
#'   `"density"` measures.
#' @param measure Character string specifying the desired measure. One of:
#'   `"survival"`, `"risk"`, `"cumulative_hazard"`, `"hazard"`, or `"density"`.
#'
#' @return Numeric value or vector of the requested measure.
#'
#' @export
convert_survival_measures <- function(survival, hazard = NULL, measure) {
  if (measure == "survival") {
    survival
  } else if (
    measure %in% c("risk", "cdf", "cumulative_distribution_function")
  ) {
    1 - survival
  } else if (measure %in% c("cumulative_hazard", "chazard")) {
    -log(survival)
  } else if (measure == "hazard") {
    hazard
  } else if (measure %in% c("density", "pdf")) {
    hazard * survival
  } else {
    cli::cli_abort(
      c(
        "The measure {.val {measure}} is not supported.",
        "i" = "Choose one of: {.val survival}, {.val risk},
               {.val cumulative_hazard}, {.val hazard}, {.val density}."
      )
    )
  }
}

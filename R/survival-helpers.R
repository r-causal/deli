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
#' @returns Numeric value or vector of the requested measure.
#'
#' @examples
#' # Risk is the complement of survival
#' convert_survival_measures(0.8, measure = "risk")
#'
#' # The density is the product of the hazard and the survival, so the
#' # "density" and "hazard" measures need the hazard as well.
#' convert_survival_measures(
#'   c(0.9, 0.8, 0.6),
#'   hazard = c(0.05, 0.06, 0.08),
#'   measure = "density"
#' )
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

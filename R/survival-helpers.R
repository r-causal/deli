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

#' The unit-time grid a supplied time design matrix is defined on
#'
#' A time design matrix models time parametrically over the unit-time intervals
#' from one to the maximum observed time, one row per interval, so the grid is
#' `seq_len(max(time))`. The two arguments that describe it may only agree with
#' it: `nrow(S)` counts its steps and `unique_times`, when supplied, names them.
#' Neither replaces it, because the same grid bins the person-periods
#' [ee_plogit()] solves on, so a fit or a prediction on another grid would be of
#' another model. [ee_plogit()] and [plogit_predict()] both come here, which is
#' what makes them agree on the grid and on the wording that refuses one.
#'
#' A maximum observed time falling between two whole times names no further whole
#' interval, and `seq_len()` truncates to the last one. The row count is compared
#' against the truncated grid, so an `S` sized to the maximum rounded up is
#' refused rather than recycled down the matrix it is added to.
#'
#' `unique_times` is compared with `all.equal()` rather than `identical()`
#' because `seq_len()` returns integers while a caller writing the grid out
#' reaches for doubles, and the two name the same grid.
#'
#' @param time The observed times, already coerced to numeric.
#' @param unique_times The `unique_times` supplied, or `NULL`.
#' @param s_rows The number of rows of the supplied time design matrix.
#' @param call The frame to report an error against, which is the function the
#'   user called.
#' @returns The integer vector `seq_len(max(time))`.
#' @noRd
plogit_unit_time_grid <- function(
  time,
  unique_times,
  s_rows,
  call = rlang::caller_env()
) {
  grid <- seq_len(max(time))

  if (length(grid) != s_rows) {
    cli::cli_abort(
      c(
        "Dimension mismatch between time intervals and {.arg S}.",
        "x" = "Found {length(grid)} unit-time intervals but {.arg S} has
               {s_rows} rows.",
        "i" = "These values must match."
      ),
      call = call
    )
  }

  agrees <- is.null(unique_times) ||
    isTRUE(all.equal(as.numeric(unique_times), as.numeric(grid)))
  if (!agrees) {
    cli::cli_abort(
      c(
        "{.arg unique_times} does not name the time grid {.arg S} is defined
         on.",
        "x" = "With {.arg S} supplied, the grid is the {length(grid)} unit-time
               interval{?s} from 1 to {length(grid)}.",
        "i" = "That grid also bins the person-periods the model is fitted on, so
               it cannot be replaced. Omit {.arg unique_times} to take it."
      ),
      call = call
    )
  }

  grid
}

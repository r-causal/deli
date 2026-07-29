#' Default names for a parameter vector
#'
#' Supplies the labels the estimators put on estimates, variance matrices, and
#' every accessor built from them. Parameters the caller named keep their names;
#' the rest are numbered by position as `theta_1` through `theta_p`.
#'
#' A partially named parameter vector is filled rather than left alone. `init =
#' c(a = 1, 2, 3)` carries the names `c("a", "", "")`, and an empty string is
#' not a usable label: it prints as a blank row and cannot be indexed by name.
#' Each empty or missing entry therefore takes the default for its position,
#' giving `c("a", "theta_2", "theta_3")`. Without that, naming one parameter
#' costs the labels of all the others.
#'
#' The result always has `p` entries, whatever the length of `nm`, so a caller
#' can label `p` parameters with it unconditionally.
#'
#' @param nm A character vector of names, or `NULL` when the parameters were
#'   given none.
#' @param p The number of parameters.
#' @returns A character vector of length `p`.
#' @noRd
default_param_names <- function(nm, p) {
  # sprintf() rather than paste0(), which treats a zero-length argument as ""
  # and would answer a zero-parameter fit with the single name "theta_".
  defaults <- sprintf("theta_%d", seq_len(p))
  if (is.null(nm)) {
    return(defaults)
  }
  # Indexing by seq_len(p) pads a short vector with NA and drops any excess, so
  # what follows works on exactly p entries.
  nm <- as.character(nm)[seq_len(p)]
  named <- !is.na(nm) & nzchar(nm)
  defaults[named] <- nm[named]
  defaults
}

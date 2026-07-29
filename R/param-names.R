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

#' Parameter names for a fitted estimator
#'
#' Resolves the two channels a caller has for labeling parameters, in order.
#' Names on `init` come first. Where `init` carries none at all, the row names
#' of the estimating functions are read instead, which is the only per-row
#' channel a `stacked_equations` closure has: the estimator sees an arbitrary
#' function returning a matrix, so nothing else it returns says which row is
#' which parameter.
#'
#' The two channels are not merged. A partially named `init` still has its gaps
#' filled positionally by `default_param_names()`, because naming one entry of a
#' vector the caller wrote out is a deliberate act. A partial row-name vector is
#' discarded whole instead, and reaches `default_param_names()` as `NULL`, so
#' the positional fill never sees one. The difference is where the partial
#' vector comes from: `rbind()` pads the rows of an unlabeled block with empty
#' strings, and `t(X * resid)` on a design whose intercept column has no name
#' produces the same shape, so an incomplete set of row names is usually an
#' accident of how the stack was built rather than a statement about any one
#' parameter. Filling it would attach a stray label or two to an otherwise
#' numbered fit.
#'
#' @param init_names The names of the estimator's `init`, or `NULL`.
#' @param evald The estimating functions evaluated at the solved values.
#' @param p The number of parameters.
#' @returns A character vector of length `p`.
#' @noRd
resolve_param_names <- function(init_names, evald, p) {
  default_param_names(init_names %||% psi_param_names(evald, p), p)
}

#' Read parameter names off an estimating-function evaluation
#'
#' Returns the row names only when they label every parameter distinctly: one
#' name per parameter, none empty, none missing, and no two alike. Anything else
#' is `NULL`. The length test is what keeps an over-identified GMM system out,
#' where the rows are moment conditions and outnumber the parameters, so its
#' labels describe the equations rather than the estimates.
#'
#' Duplicates are discarded for the same reason an incomplete set is: this
#' channel carries labels the caller never typed. `rbind()` names a row after
#' the variable that supplied it, so stacking two blocks that each begin with a
#' variable of the same name yields a repeated label from a stack that looks
#' unremarkable. A repeated label is worse than no label, because `coef()` then
#' shows one name twice and `confint(m)["mu", ]` silently returns whichever row
#' comes first. Names on `init` are left alone by this rule: a caller who writes
#' `c(mu = 0, mu = 1)` has typed the repetition out.
#'
#' @param evald The estimating functions evaluated at the solved values.
#' @param p The number of parameters.
#' @returns A character vector of length `p`, or `NULL`.
#' @noRd
psi_param_names <- function(evald, p) {
  nm <- rownames(evald)
  if (is.null(nm) || length(nm) != p) {
    return(NULL)
  }
  nm <- as.character(nm)
  if (anyNA(nm) || !all(nzchar(nm)) || anyDuplicated(nm) > 0L) {
    return(NULL)
  }
  nm
}

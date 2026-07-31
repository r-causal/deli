# ---- parameter-name conditions -----------------------------------------------
# The abort raised while labeling a parameter vector carries a condition class,
# so a caller or a test can match on the class rather than on the prose, as the
# initial-value failures in `R/validation.R` do.
#
#   deli_param_name_collision
#     The positional fill would give an unnamed parameter a name another
#     parameter in the same vector already carries. Raised by
#     default_param_names(), which is the one place a label no caller typed is
#     created, and so the one place such a repetition can come from.

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
#' What the fill will not do is synthesize a name another parameter already
#' carries. `init = c(theta_2 = 0, 1)` carries the names `c("theta_2", "")`, and
#' the default for the empty entry is `theta_2`, the name position 1 was given.
#' A repeated label is worse than no label, because `coef()` then shows one name
#' twice and `confint(m)["theta_2", ]` silently returns whichever row comes
#' first, so that fill is refused rather than performed. Only a repetition the
#' fill would create is refused: one the caller typed out stands, for the reason
#' `psi_param_names()` gives below. A name matching the default for its own
#' position repeats nothing and is filled around as usual.
#'
#' The result always has `p` entries, whatever the length of `nm`, so a caller
#' can label `p` parameters with it unconditionally.
#'
#' @param nm A character vector of names, or `NULL` when the parameters were
#'   given none.
#' @param p The number of parameters.
#' @param error_call The frame to report the refusal against. This function is
#'   several frames below whatever the caller reached and names its own
#'   parameters, so a fit passes the frame of the method that was called.
#'   `NULL`, the default, reports no call, which is what the print and summary
#'   paths leave: they reach the fill from a parameter vector renamed after the
#'   fit, and the frame that displays a name is not the one that set it.
#' @returns A character vector of length `p`. Raises an error carrying the class
#'   `deli_param_name_collision` where the fill would repeat a name the vector
#'   already holds.
#' @noRd
default_param_names <- function(nm, p, error_call = NULL) {
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
  # Judged over the whole vector before any of it is filled, so the report
  # covers every position the numbering would repeat at rather than however many
  # it had reached, and so a refused vector is left exactly as it arrived.
  repeated <- which(!named & defaults %in% nm[named])
  if (length(repeated) > 0L) {
    # The positions are written out as strings because cli reads a single number
    # as the quantity to pluralize on rather than counting the values it was
    # given, so one collision at position 2 would be described in the plural.
    positions <- as.character(repeated)
    cli::cli_abort(
      c(
        "The unnamed parameters cannot be numbered without repeating a name.",
        "x" = "The numbering would put {.val {defaults[repeated]}} at
               position{?s} {positions}, and {?that name is/those names are}
               already in use.",
        "i" = "Name every parameter, or use names outside the
               {.code theta_<i>} pattern the numbering follows.",
        "i" = "Names inherited from a previous fit's estimates collide this
               way; drop them with {.fn unname}, or name the remaining
               parameters."
      ),
      class = "deli_param_name_collision",
      call = error_call
    )
  }
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
#' A formula fit reaches this with names on `init` more often than the caller
#' typed any: `prepare_formula_psi()` fills an unnamed `init` from the model
#' matrix through `design_param_names()` below, since the column headings it
#' knows are the ones no estimating equation is given.
#'
#' @param init_names The names of the estimator's `init`, or `NULL`.
#' @param evald The estimating functions evaluated at the solved values.
#' @param p The number of parameters.
#' @param error_call The frame to report a refused fill against, passed through
#'   to `default_param_names()`.
#' @returns A character vector of length `p`.
#' @noRd
resolve_param_names <- function(init_names, evald, p, error_call = NULL) {
  default_param_names(
    init_names %||% psi_param_names(evald, p),
    p,
    error_call = error_call
  )
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

#' Positional names for one block of coefficients
#'
#' Labels the `p` coefficients of a regression block `prefix_1` through
#' `prefix_p`. The built-in estimating equations use it for the nuisance blocks
#' they stack beneath the parameters a caller came for, taking the prefix from
#' the design argument that supplied the block, so `W_2` names the coefficient
#' on the second column of `W`.
#'
#' Positional is all that is available and all that is claimed. Every design
#' argument runs through `coerce_design()`, which drops its column headings, so
#' an estimating equation knows how many coefficients a block holds and which
#' argument they belong to but not what the caller called them.
#'
#' Labeling these blocks at all is what lets the interesting parameters keep
#' their names. `psi_param_names()` above discards a row-name vector that does
#' not label every parameter, so leaving a nuisance block blank would cost the
#' causal estimates beside it their labels as well.
#'
#' Which built-in equations name their rows follows from that. An equation names
#' them when at least one row carries a label a design position cannot supply:
#' a quantity the equation defines, such as `ACE` or `lambda` or
#' `log_inv_scale`, or a symbol the literature gives a whole block, such as
#' `SNM phi` or `MSM alpha`. The design blocks stacked beside that row are then
#' named here so the set is complete. An equation whose parameters are all
#' coefficients on a design has no such row, and positional labels alone would
#' say no more than the numbering they replace, so `ee_regression()`,
#' `ee_plogit()`, `ee_glm()` under every distribution but gamma and negative
#' binomial, and `ee_aft()` under the exponential distribution leave their rows
#' unnamed. So do the single-quantity helpers built to be stacked more than
#' once, `ee_mean()` and `ee_percentile()` and the two `_ed()` equations among
#' them, because a fixed label repeated across a stack costs that stack its
#' names.
#'
#' Four more returns follow from the same rule. `ee_glm()` under the gamma and
#' negative binomial distributions, `ee_tobit()`, and `ee_beta_regression()`
#' each hold one row more than their design has columns, and that row is a
#' parameter of the outcome distribution rather than a coefficient on a design
#' column, which is the property that earns `ee_aft()`'s `log_inv_scale` its
#' label. Each names that row `log_shape`, `log_dispersion`, `log_sigma`, or
#' `log_phi`, with the design block named here beside it.
#'
#' @param prefix The name of the block, ordinarily the design argument that
#'   supplied it.
#' @param p The number of coefficients in the block.
#' @returns A character vector of length `p`.
#' @noRd
block_param_names <- function(prefix, p) {
  sprintf("%s_%d", prefix, seq_len(p))
}

#' The parameter an estimating equation appends past its design coefficients
#'
#' Answers, for the built-in equations the formula interface can drive, whether
#' the equation estimates one parameter beyond a coefficient per design column
#' and what that parameter is called. `NULL` for every equation whose parameters
#' are all design coefficients, and for every equation not listed.
#'
#' The table is written out rather than derived from what the equations return.
#' Two things it has to know cannot be read off a return. The first is where in
#' `theta` an appended parameter sits: `ee_gformula()` labels its rows
#' `c("causal_mean", "X_1", ...)` for one parameter more than the design has
#' columns, the same count as `ee_tobit()`, with the extra parameter at the
#' front rather than the back, so the row names alone do not tell the two
#' layouts apart. The second is that the answer is needed before the equation
#' has evaluated successfully at all, since the diagnostic in
#' `abort_formula_auto_init()` reports on an automatic `init` that the equation
#' rejected.
#'
#' The names here must match the ones those equations write on their own rows,
#' and `test-param-names-builtin.R` pins both ends.
#'
#' @param .ee The estimating-equation function passed to the formula interface.
#' @param ee_args The evaluated `...` arguments forwarded to it.
#' @returns A string, or `NULL`.
#' @noRd
appended_param_name <- function(.ee, ee_args) {
  if (!is.function(.ee)) {
    return(NULL)
  }
  dist <- ee_args[["distribution"]]
  if (!is.character(dist) || length(dist) != 1L) {
    dist <- ""
  }
  dist <- tolower(dist)

  if (identical(.ee, ee_glm)) {
    if (dist == "gamma") {
      return("log_shape")
    }
    if (dist %in% c("negative_binomial", "nb")) {
      return("log_dispersion")
    }
    # Every other distribution, the tweedie among them, estimates one
    # coefficient per design column.
    return(NULL)
  }
  if (identical(.ee, ee_aft)) {
    # Sigma is fixed at 1 under the exponential, so that distribution alone
    # appends nothing. An unreadable distribution is treated the same way,
    # since a wrong label is worse than the labels the equation writes on its
    # own rows, and worse than numbering where it writes none.
    if (dist %in% c("", "exponential")) {
      return(NULL)
    }
    return("log_inv_scale")
  }
  if (identical(.ee, ee_tobit)) {
    return("log_sigma")
  }
  if (identical(.ee, ee_beta_regression)) {
    return("log_phi")
  }
  NULL
}

#' Names for an unnamed `init` on the formula interface
#'
#' The formula interface knows one thing no estimating equation does: what the
#' caller called the columns of the design. `coerce_design()` drops the column
#' headings on the way in, so an equation can label a design block only by
#' position and the row-name channel cannot recover them. A formula fit whose
#' `init` carries no names of its own therefore takes them from the model
#' matrix, and from the parameter the equation appends beyond it.
#'
#' Only two lengths are named: one parameter per design column, and one per
#' design column plus the appended parameter. Any other length describes a
#' parameter vector this function cannot account for, and naming it would
#' attach the intercept's label to whatever the first parameter happens to be,
#' so it is left alone. Leaving it alone hands the labeling to the row names of
#' the estimating functions, which is how a formula fit of `ee_gformula()` comes
#' back labeled even though its extra parameter leads rather than trails. The
#' parameters are numbered only where those row names do not label every one.
#'
#' @param design_names The column names of the model matrix.
#' @param appended The parameter the estimating equation appends, or `NULL`.
#' @param p The number of parameters.
#' @returns A character vector of length `p`, or `NULL`.
#' @noRd
design_param_names <- function(design_names, appended, p) {
  if (is.null(design_names)) {
    return(NULL)
  }
  k <- length(design_names)
  if (p == k) {
    return(design_names)
  }
  if (!is.null(appended) && p == k + 1L) {
    return(c(design_names, appended))
  }
  NULL
}

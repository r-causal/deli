#' Aggregate estimating function contributions by group
#'
#' @description
#' Collapses unit-level estimating function contributions into group-level
#' contributions for clustered or grouped data. Uses an independent working
#' correlation structure (summing within groups).
#'
#' This function mirrors `aggregate_efuncs()` in Python delicatessen, so code
#' translated from Python can keep its shape. There is no base R equivalent that
#' operates on estimating function contributions, so this is the interface for
#' them in deli as well.
#'
#' @details
#' This function should be called inside the `psi` function after computing
#' unit-level estimating equations but before returning them to
#' [MEstimator()]. This changes the effective sample size used by the
#' empirical sandwich variance estimator.
#'
#' @param est_funcs A p-by-n matrix of estimating function contributions,
#'   where p is the number of parameters and n is the number of observations.
#'   A length-n vector is treated as a single parameter observed across n
#'   observations, matching a 1-by-n matrix.
#' @param group A vector of length n identifying the group (cluster) for
#'   each observation.
#'
#' @returns A p-by-m matrix, where m is the number of unique groups. Row names
#'   are those of `est_funcs`, since the rows are the same parameters. Columns
#'   are ordered by the sorted unique values of `group` and are labeled with
#'   those values, as character. A factor `group` is coerced with `as.vector()`
#'   to its character labels before sorting, so its columns sort lexically by
#'   label rather than by factor-level order and carry those labels; a level
#'   with no observations contributes no column and so no label.
#'
#' @examples
#' # Fifty clusters of four observations, sharing a cluster-level shift in y
#' set.seed(42)
#' n <- 200
#' group <- rep(1:50, each = 4)
#' cluster_effect <- rnorm(50, sd = 2)
#' y <- cluster_effect[group] + rnorm(n)
#'
#' psi <- function(theta) aggregate_efuncs(ee_mean(theta, y = y), group = group)
#'
#' m <- m_estimate(stacked_equations = psi, init = mean(y))
#'
#' # Cluster-robust standard error, larger than the naive independence version
#' sqrt(diag(vcov(m)))
#'
#' @export
aggregate_efuncs <- function(est_funcs, group) {
  # Summing within a group is a linear operation, so a derivative aggregates the
  # way the values do: the derivative of a group's summed contribution is the
  # sum of the derivatives of the contributions in it. Aggregating each slot of
  # a tangent-carrying return is therefore exact rather than an approximation of
  # it. The recursion has to intercept ahead of everything below, which asks a
  # plain matrix to hold a derivative it has nowhere to put.
  if (is_tangent_container(est_funcs)) {
    parts <- pt_arrays(est_funcs)
    # The two slots are aggregated separately, so a tangent that broadcasts over
    # the primal has to reach the primal's length and shape first.
    tangent <- pt_recycle_tangent(parts$primal, parts$tangent)
    return(primal_tangent_array(
      aggregate_efuncs(parts$primal, group),
      aggregate_efuncs(tangent, group)
    ))
  }

  # A 1-D input is a single parameter observed across n units, matching the
  # Python reference, so it becomes a 1-by-n row rather than the n-by-1 column
  # that `as.matrix` would otherwise produce from a vector.
  if (is.null(dim(est_funcs))) {
    est_funcs <- matrix(est_funcs, nrow = 1)
  } else {
    est_funcs <- as.matrix(est_funcs)
  }
  group <- as.vector(group)

  # Get dimensions
  n_prm <- nrow(est_funcs)
  n_obs <- ncol(est_funcs)

  # Validate input dimensions
  if (length(group) != n_obs) {
    cli::cli_abort(
      c(
        "Length of {.arg group} must match the number of columns in
         {.arg est_funcs}.",
        "x" = "{.arg group} has {length(group)} element{?s} but
               {.arg est_funcs} has {n_obs} column{?s}."
      )
    )
  }

  # Map groups to compact integer indices. Sort the unique groups so that the
  # aggregated columns come back in sorted unique-group order, matching the
  # Python reference (which relies on np.unique).
  unique_groups <- sort(unique(group))
  m <- length(unique_groups)
  group_idx <- match(group, unique_groups)

  # Aggregate by summing within groups using rowsum on transposed matrix
  # t(est_funcs) is n-by-p, rowsum sums rows by group -> m-by-p
  aggregated <- rowsum(t(est_funcs), group_idx, reorder = TRUE)

  # Return as p-by-m matrix. The rows are the same parameters as the rows of
  # the input, so any labels on them survive. The columns are the groups, and
  # they are labeled with the group values themselves. What rowsum() labels its
  # own rows with is the compact indices it was handed, which stand for nothing
  # outside this function, so those labels are replaced rather than passed on.
  out <- unname(t(aggregated))
  rownames(out) <- rownames(est_funcs)
  colnames(out) <- as.character(unique_groups)
  out
}

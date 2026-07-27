#' Aggregate estimating function contributions by group
#'
#' Collapses unit-level estimating function contributions into group-level
#' contributions for clustered or grouped data. Uses an independent working
#' correlation structure (summing within groups).
#'
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
#' @returns A p-by-m matrix, where m is the number of unique groups. Columns
#'   are ordered by the sorted unique values of `group`. A factor `group` is
#'   coerced with `as.vector()` to its character labels before sorting, so its
#'   columns sort lexically by label rather than by factor-level order.
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
#' m <- MEstimator(stacked_equations = psi, init = mean(y)) |>
#'   estimate()
#'
#' # Cluster-robust standard error, larger than the naive independence version
#' sqrt(diag(vcov(m)))
#'
#' @export
aggregate_efuncs <- function(est_funcs, group) {
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

  # Return as p-by-m matrix (strip names from rowsum)
  unname(t(aggregated))
}

#' GMM Estimator
#'
#' S7 class for Generalized Method of Moments (GMM) estimation via
#' minimization of estimating equations with empirical sandwich variance
#' estimation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of estimating equations and n is the number of observations.
#'   The number of equations p must be greater than or equal to the number of
#'   parameters (length of `init`). Row names on that matrix name the parameters
#'   when `init` has none and there is exactly one row per parameter, so an
#'   over-identified system is numbered instead; see [estimate()].
#' @param init Numeric vector of initial parameter values for the minimization
#'   algorithm. Names on it label the parameters and take precedence over the
#'   row names of `stacked_equations`.
#' @inheritParams gmm_estimate
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction.
#' @param overid_maxiter Integer maximum iterations for the two-step iterative
#'   procedure for over-identified problems. Default `200L`. The update converges
#'   linearly rather than quadratically, so a well-identified system commonly
#'   needs tens of passes to reach `overid_tolerance` and a weakly identified one
#'   can need hundreds.
#' @param overid_tolerance Numeric tolerance for convergence of the two-step
#'   iterative procedure for over-identified problems. Default `1e-9`.
#' @param summed_equations A function that takes a numeric vector `theta` and
#'   returns the length-p vector of row sums of `stacked_equations` at `theta`,
#'   or `NULL` (default) to derive those sums from the full p-by-n return.
#'
#'   Two of the three things a GMM fit does with the estimating functions want
#'   nothing but those sums. The objective is a quadratic form in the mean
#'   moments, \eqn{\bar{g}(\theta)' W \bar{g}(\theta)}, so every evaluation the
#'   minimizer makes reduces the whole p-by-n matrix to the reduction's own
#'   output; the bread is the Jacobian of the same sums. Both take a supplied
#'   reduction and skip the matrix.
#'
#'   The third does not. The two-step weight matrix is the inverse of the
#'   covariance of the moment conditions, which is a cross-product of the
#'   per-observation contributions, so each pass of the update over an
#'   over-identified system evaluates `stacked_equations` in full whatever this
#'   property holds. So do the validation at the starting values and the meat.
#'   A just-identified fit runs no weight update, so it makes exactly those two
#'   full evaluations.
#'
#'   Under `deriv_method = "exact"` the reduction is called with a
#'   tangent-carrying `theta`, so it must be written in operations that carry
#'   derivatives: `t(X) %*% r` does, and [base::crossprod()] does not. See
#'   [auto_differentiation()] for which operations carry a tangent and where.
#'
#'   Anything that is neither `NULL` nor a function is refused here, with an
#'   error carrying the class `deli_summed_equations_error`. So is a return at
#'   the estimated values that is not numeric or holds fewer values than the
#'   system has moment conditions, which [estimate()] reads.
#'
#'   [gmm_estimate()] and [m_estimate()] do not offer it, for the reason
#'   [MEstimator()] gives: both build `stacked_equations` themselves, so the
#'   moment conditions a reduction would have to match are ones the caller never
#'   writes.
#' @param check_summed_equations Logical. When `TRUE` (default) and
#'   `summed_equations` was supplied, its value at the estimated values is
#'   compared against the row sums of the one full evaluation the meat is built
#'   from, and a disagreement raises an error carrying the class
#'   `deli_summed_equations_disagree`. See [MEstimator()] for what the
#'   comparison catches, what it cannot, and what setting it to `FALSE` leaves
#'   behind.
#'
#' @returns A `GMMEstimator` S7 object. Call [estimate()] to minimize the
#'   estimating equations and compute the sandwich variance.
#'
#' @section Moment quality of an over-identified fit:
#' A just-identified system has as many moment conditions as parameters, so the
#' moments vanish at a solution and the size of what is left over says whether the
#' fit succeeded. An over-identified system has no such reading: no value of the
#' parameters drives every condition to zero, and a residual moment is expected
#' rather than diagnostic. Hansen's J-statistic is the reading that is available
#' there. It is n times the GMM objective at the minimum,
#' \eqn{J = n \bar{g}(\hat{\theta})' W \bar{g}(\hat{\theta})}, where
#' \eqn{\bar{g}} averages the moment conditions over the observations and
#' \eqn{W} is the weight matrix the fit finished with. Under correct
#' specification it is asymptotically chi-squared on as many degrees of freedom as
#' the system has moment conditions beyond parameters, so its size can be judged
#' against a reference distribution rather than against the scale of the data.
#'
#' [estimate()] records it in the `j_statistic` property of an over-identified
#' fit, and [`summary()`][deli-display] reports it with its degrees of freedom and
#' its P-value. A just-identified fit has no degrees of freedom left over and
#' leaves the property `NULL`; its moments are judged directly instead, as
#' [estimate()] describes. A `subset` fit holds the parameters outside the subset
#' at their initial values rather than estimating them, which the reference
#' distribution does not allow for, so it is left `NULL` too.
#'
#' A P-value the reference distribution all but rules out warns with the class
#' `deli_gmm_moments_rejected`, which usually means the moment conditions cannot
#' all hold at one value of the parameters. The weight matrix is what makes J
#' comparable across problems, so the warning is raised only where the two-step
#' update settled: a fit that exhausted `overid_maxiter` has already warned about
#' that, and its J has no reference distribution to be judged against. The
#' property still records the statistic in that case, as it does for
#' `overid_maxiter = 0`, which leaves the identity weight matrix in place and so
#' leaves J an unstandardized sum of squared moments.
#'
#' The reading J cannot make is the opposite failure. Moment conditions that are
#' linearly dependent, one of them repeating what the others already say, leave
#' the covariance the weight matrix inverts singular, and the update falls
#' through to the pseudo-inverse; the fit that comes back is the fit of the
#' independent conditions alone. J is silent about it, because a condition the
#' others account for agrees with them wherever the parameters sit and so adds
#' nothing for J to measure, which drives J toward zero rather than away from
#' it. That case warns with the class `deli_gmm_moments_dependent` instead,
#' naming the conditions the factorization found redundant.
#'
#' @export
#' @examples
#' # The constructor builds the estimator and `estimate()` solves it, so an
#' # object that has not been through `estimate()` reports only what it was
#' # given. `gmm_estimate()` does both steps in one call.
#' y <- c(1, 2, 3, 4, 5)
#' psi <- function(theta) {
#'   matrix(y - theta[1], nrow = 1)
#' }
#' GMMEstimator(stacked_equations = psi, init = 0)
#'
#' # One moment condition for one parameter is just-identified, so the minimizer
#' # lands where `MEstimator()` would have found the root.
#' GMMEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#'
#' # A Poisson mean is identified twice over, by the mean and by the variance,
#' # so these two moment conditions estimate one parameter and the system is
#' # over-identified. That is the case `MEstimator()` cannot solve, and the case
#' # the `overid_maxiter` and `overid_tolerance` properties govern: they stop
#' # the two-step weight matrix update that reconciles the two conditions.
#' set.seed(42)
#' counts <- rpois(200, lambda = 3)
#'
#' psi_pois <- function(theta) {
#'   rbind(
#'     counts - theta[1],
#'     (counts - theta[1])^2 - theta[1]
#'   )
#' }
#'
#' g <- GMMEstimator(stacked_equations = psi_pois, init = 1) |>
#'   estimate()
#'
#' # With more conditions than parameters neither is solved exactly. The weight
#' # matrix is what decides how the disagreement between them is split.
#' summary(g)
GMMEstimator <- new_class(
  "GMMEstimator",
  parent = deli_estimator,
  properties = list(
    overid_maxiter = new_property(
      class = class_integer,
      default = 200L
    ),
    overid_tolerance = new_property(
      class = class_numeric,
      default = 1e-9
    ),
    weight_matrix = new_property(
      class = NULL | class_double,
      default = NULL
    ),
    j_statistic = new_property(
      class = NULL | class_double,
      default = NULL
    )
  ),
  constructor = function(
    stacked_equations,
    init,
    subset = NULL,
    finite_correction = NULL,
    overid_maxiter = 200L,
    overid_tolerance = 1e-9,
    summed_equations = NULL,
    check_summed_equations = TRUE
  ) {
    check_estimator_init(init)
    check_finite_correction(finite_correction)
    check_estimator_subset(subset, length(init))
    check_overid_scalar(overid_maxiter, "overid_maxiter")
    check_overid_scalar(overid_tolerance, "overid_tolerance")
    check_summed_function(summed_equations)
    check_summed_switch(check_summed_equations)

    init_names <- names(init)
    init <- as.numeric(init)
    names(init) <- init_names

    # Sort subset if provided
    if (!is.null(subset)) {
      subset <- sort(as.integer(subset))
    }

    new_object(
      S7_object(),
      stacked_equations = stacked_equations,
      summed_equations = summed_equations,
      check_summed_equations = check_summed_equations,
      init = init,
      subset = subset,
      finite_correction = finite_correction,
      overid_maxiter = as.integer(overid_maxiter),
      overid_tolerance = overid_tolerance,
      n_params = length(init)
    )
  }
)

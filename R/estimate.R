#' Estimate parameters and sandwich variance
#'
#' Solves the estimating equations for the parameter vector `theta` and
#' computes the empirical sandwich variance estimator.
#'
#' @details
#' The estimates and every matrix built from them carry parameter names, taken
#' from the first of two channels that supplies them. Names on `init` come
#' first: a parameter the caller named keeps that name, and one left unnamed is
#' numbered by position, so `init = c(mu = 0, 1)` gives `c("mu", "theta_2")`.
#' Where `init` carries no names at all, the row names of the estimating
#' functions are read instead. That is how a `stacked_equations` function names
#' the parameters it defines, since the estimator otherwise sees only an opaque
#' closure returning a matrix. The formula interface names `init` from the model
#' matrix columns where it can account for its length, one parameter per design
#' column or one more for the parameter the equation appends, so the row names
#' carry the function interface and every formula fit of another length.
#' [ee_gformula()] is one of those: its extra parameter leads rather than
#' trails, so the formula interface leaves its `init` unnamed and the fit takes
#' the labels the equation writes on its own rows.
#'
#' Row names are read only when they label every parameter distinctly: one name
#' per parameter, none empty, none missing, and no two alike. An incomplete or
#' repetitive set is discarded rather than patched up, because it is usually an
#' accident of how the stack was built. `rbind()` pads the rows of an unlabeled
#' block with empty strings, and `t(X * resid)` on a design whose intercept
#' column has no name produces the same shape; `rbind()` also names each row
#' after the variable that supplied it, so two blocks that each begin with a
#' variable of the same name repeat a label. The count has to match as well,
#' which is what keeps an over-identified `GMMEstimator` out: its rows are
#' moment conditions and outnumber the parameters, so their labels describe the
#' equations. Where neither channel applies, the parameters are numbered
#' `theta_1` through `theta_p`.
#'
#' Row names survive exact differentiation. Assigning them is ignored while a
#' value carries derivatives, since the labels are read from the plain
#' evaluation at the solved values. `rownames<-` hands its work to `dimnames<-`,
#' which is a generic, so deli registers the setter there and an estimating
#' function can label its rows from anywhere.
#'
#' Many of the built-in estimating equations name their own rows, so a fit that
#' passes one of them a wholly unnamed `init` comes back labeled rather than
#' numbered. Each documents its labels under **Value**. The rule above still
#' governs a stack built from them: two blocks that name the same parameter, as
#' two [ee_ipw()] blocks do, repeat a label and the fit is numbered instead, and
#' one named block stacked with an unnamed one is incomplete and is numbered as
#' well. Name `init` where a stack needs labels the blocks cannot agree on.
#'
#' @param object An `MEstimator` or `GMMEstimator` object.
#' @param solver Character string specifying the solver algorithm, or a custom
#'   function. When `NULL` (default), uses `"rootSolve"` for `MEstimator`
#'   ([rootSolve::multiroot()]) and `"BFGS"` for `GMMEstimator`
#'   ([stats::optim()]). Other options for `MEstimator`: `"lm"`, the
#'   Levenberg-Marquardt algorithm ([minpack.lm::nls.lm()]), which mirrors the
#'   default solver of Python `delicatessen`
#'   (`scipy.optimize.root(method = "lm")`); and `"nleqslv"` (uses
#'   [nleqslv::nleqslv()]). A custom function must accept
#'   `stacked_equations` and `init` arguments and return the solved theta
#'   vector. Whichever solver is used, the point it returns is judged against
#'   the estimating equations themselves and a warning is raised when they are
#'   not solved there. A custom function reports no status of its own, so its
#'   point is judged exactly as a built-in solver's is. Two `GMMEstimator`
#'   fits are exceptions. An over-identified one cannot drive every moment to
#'   zero, so its moments are read against the J-statistic instead, which
#'   [GMMEstimator()] describes. A `subset` one holds some parameters fixed while
#'   still summing every equation into the objective, so it is judged neither way
#'   and does not warn; inspect `rowSums()` of the estimating functions at the
#'   returned values. [rootSolve::multiroot()] cannot run inside itself,
#'   so a fit whose estimating function fits a second M-estimator cannot leave
#'   both on the default solver; that is refused rather than attempted, and
#'   naming `"nleqslv"` for either of the two fits resolves it.
#' @param maxiter Integer maximum iterations for the solver (default 5000).
#' @param tolerance Numeric tolerance for the solver (default 1e-9).
#' @param deriv_method Character string for the derivative method used to
#'   compute the bread matrix. One of `"capprox"` (central difference),
#'   `"fapprox"` (forward difference), `"bapprox"` (backward difference), or
#'   `"exact"` (forward-mode automatic differentiation). Default `"capprox"`.
#'   Exact differentiation removes the finite-difference step size but requires
#'   that the estimating equation is composed of operations the autodiff
#'   supports (see [auto_differentiation()]). For post-estimation transforms and
#'   the survival prediction helpers, exact differentiation is also available
#'   through [delta_method()].
#' @param dx Numeric step size for numerical differentiation; ignored when
#'   `deriv_method = "exact"` (default 1e-9). Must be a single positive finite
#'   number, which is checked whichever `deriv_method` is in force. The step is
#'   absolute and is floored at the floating-point resolution of each estimate,
#'   so a large parameter magnitude cannot silently reduce it to nothing; see
#'   [approx_differentiation()].
#' @param allow_pinv Logical. Use pseudo-inverse if bread is singular? Default
#'   `TRUE`.
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A modified `MEstimator` object with populated `theta`, `bread`,
#'   `meat`, `variance`, and `asymptotic_variance` properties.
#'
#' @export
#' @examples
#' psi <- function(theta) {
#'   y <- c(1, 2, 3, 4, 5)
#'   matrix(y - theta[1], nrow = 1)
#' }
#' m <- MEstimator(stacked_equations = psi, init = 0) |>
#'   estimate()
#' coef(m)
estimate <- new_generic(
  "estimate",
  "object",
  function(
    object,
    solver = NULL,
    maxiter = 5000,
    tolerance = 1e-9,
    deriv_method = "capprox",
    dx = 1e-9,
    allow_pinv = TRUE,
    ...
  ) {
    S7::S7_dispatch()
  }
)

method(estimate, MEstimator) <- function(
  object,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_dx(dx)
  # One fit evaluates the estimating function many times, so an estimating
  # function that warns raises the same warning repeatedly for one operation.
  # See R/conditions.R. The scope wraps a worker rather than this body because a
  # calling handler cannot be installed for a frame that is already running, and
  # check_dots_empty() has to stay here, where caller_env() is the user's frame.
  without_repeated_warnings(estimate_m_estimator(
    object,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  ))
}

#' Solve an `MEstimator` and assemble its sandwich variance
#'
#' The body of the [estimate()] method for [MEstimator()], separated from it so
#' that the method can wrap the whole of it in `without_repeated_warnings()`.
#'
#' @inheritParams estimate
#' @returns The `MEstimator` with its estimated properties populated.
#' @noRd
estimate_m_estimator <- function(
  object,
  solver,
  maxiter,
  tolerance,
  deriv_method,
  dx,
  allow_pinv
) {
  # Default solver for MEstimator is rootSolve
  if (is.null(solver)) {
    solver <- "rootSolve"
  }
  stacked_equations <- object@stacked_equations
  init <- object@init
  subset <- object@subset
  finite_correction <- object@finite_correction

  # Validate the estimating-function return at the initial values before
  # handing off to the solver, so a malformed return produces an informative
  # error instead of an opaque failure inside the solver. The summed equations
  # there are kept for judging the returned point, which needs to know which
  # equations the solver moved; nothing else evaluates the estimating function
  # at the starting values, so keeping them costs no evaluation. What the
  # validation judges is the estimating function the caller supplied, so it
  # names the frame of the method that was called, as the solver errors below
  # do.
  init_score <- equation_scores(eval_psi_at_init(
    stacked_equations,
    init,
    error_call = rlang::caller_env()
  ))

  # Build the summed EE function for root-finding
  summed_ee <- function(theta) {
    # If subset, expand theta into full parameter vector
    if (!is.null(subset)) {
      full_theta <- init
      full_theta[subset] <- theta
    } else {
      full_theta <- theta
    }

    ef <- stacked_equations(full_theta)
    if (is.null(dim(ef))) {
      s <- sum(ef)
    } else {
      s <- rowSums(ef)
    }

    # If subset, return only the subset equations
    if (!is.null(subset)) {
      s <- s[subset]
    }
    s
  }

  # Get initial values (possibly subset)
  if (!is.null(subset)) {
    inits <- init[subset]
  } else {
    inits <- init
  }

  # Solve
  if (is.character(solver)) {
    # The errors solve_equations() raises are about what the caller asked for,
    # so they name the frame of the method that was called rather than either of
    # the internal functions in between.
    solved <- solve_equations(
      summed_ee,
      inits,
      solver,
      maxiter,
      tolerance,
      call = rlang::caller_env()
    )
  } else if (is.function(solver)) {
    par <- solver(stacked_equations = summed_ee, init = inits)
    check_solver_return(par, length(inits))
    # A custom solver reports nothing, so its point is judged like one a
    # built-in solver declared solved.
    solved <- list(par = par, solver = "custom", warned = FALSE)
  } else {
    # This is the only user-reachable condition raised directly in the worker,
    # and the worker is internal, so the error names the frame of the method
    # that called it rather than a function no user can see.
    cli::cli_abort(
      "{.arg solver} must be a string or function.",
      call = rlang::caller_env()
    )
  }
  theta_solved <- solved$par

  # Expand theta if subset was used
  if (!is.null(subset)) {
    full_theta <- init
    full_theta[subset] <- theta_solved
  } else {
    full_theta <- theta_solved
  }

  # Get n_obs from evaluating at solved theta
  evald <- stacked_equations(full_theta)
  if (is.null(dim(evald))) {
    n_obs <- length(evald)
    # Reshape to 1-by-n matrix so tcrossprod works correctly in compute_meat
    evald <- matrix(evald, nrow = 1)
  } else {
    n_obs <- ncol(evald)
  }

  # Compute sandwich components
  bread <- compute_bread(stacked_equations, full_theta, deriv_method, dx) /
    n_obs

  # Account for the solve, unless the solver reported a failure the caller has
  # already been warned about. Everything here waits for the bread: the Jacobian
  # is what measures the distance still to travel, a solver whose own convergence
  # test is relative to its starting residuals can report success without moving
  # at all, and a solver that reports failure cannot say whether the search or
  # the problem is what failed. A subset fit solves only the subset equations for
  # the subset parameters, so only that block is read.
  #
  # The point is read once and reported through exactly one of two branches, so
  # the solve raises at most one warning. Where the solver reported a failure the
  # reading decides the wording rather than whether to speak, since a point that
  # solves the equations says the search is not what went wrong.
  if (!solved$warned) {
    judged <- subset %||% seq_len(nrow(evald))
    judged_bread <- bread[judged, judged, drop = FALSE]
    unsolved <- unsolved_point(
      evald[judged, , drop = FALSE],
      full_theta[judged],
      judged_bread,
      init_score = init_score[judged],
      init = init[judged]
    )
    if (isTRUE(solved$not_converged)) {
      warn_solver_failure(solved, judged_bread, unsolved)
    } else if (!is.null(unsolved)) {
      warn_unsolved(solved, unsolved)
    }
  }

  meat_mat <- compute_meat(evald) / n_obs
  meat_mat <- finite_sample_correction(
    meat_mat,
    n_obs,
    length(full_theta),
    finite_correction
  )
  asymp_var <- build_sandwich(bread, meat_mat, allow_pinv)

  if (is.null(asymp_var)) {
    var_mat <- NULL
  } else {
    var_mat <- asymp_var / n_obs
  }

  # Apply parameter names, from `init` where it has any and otherwise from the
  # row names of the evaluation already made above. No second call to the
  # estimating functions is needed to read them. A fill this refuses is about
  # the `init` the caller wrote, so it names the frame of the method that was
  # called rather than either of the naming helpers below it.
  param_names <- resolve_param_names(
    names(object@init),
    evald,
    length(full_theta),
    error_call = rlang::caller_env()
  )
  names(full_theta) <- param_names
  dimnames(bread) <- list(param_names, param_names)
  dimnames(meat_mat) <- list(param_names, param_names)
  if (!is.null(asymp_var)) {
    dimnames(asymp_var) <- list(param_names, param_names)
  }
  if (!is.null(var_mat)) {
    dimnames(var_mat) <- list(param_names, param_names)
  }

  # Update the object
  object@theta <- full_theta
  object@n_obs <- as.integer(n_obs)
  object@bread <- bread
  object@meat <- meat_mat
  object@asymptotic_variance <- asymp_var
  object@variance <- var_mat

  object
}

# ---- solver conditions -------------------------------------------------------
# The conditions raised about a solve carry a condition class, so a caller or a
# test can match on the class rather than on the prose, as the exact-mode aborts
# in `R/autodiff.R` do.
#
#   deli_solver_not_converged
#     The solver stopped without solving the estimating equations, either
#     because it reported a failure of its own or because the returned point
#     does not solve them. Every solver branch raises this one class. An
#     M-estimation solve reports itself through exactly one of three places, so
#     it raises the class at most once: an exhausted iteration budget or a
#     solver that failed outright warns from solve_equations(), a solver whose
#     own convergence test failed warns from warn_solver_failure(), and a solver
#     that reported success has only its returned point judged. A just-identified
#     GMM solve judges its moments only where the minimizer reported success, so
#     it too raises it at most once; an over-identified solve calls the minimizer
#     repeatedly and can raise it once per pass, though the scope described in
#     R/conditions.R delivers the passes that report the same thing as one
#     warning.
#
#   deli_gmm_moments_rejected
#     The Hansen J-statistic of an over-identified GMM fit is so large that the
#     chi-squared reference distribution all but rules the fit out, which usually
#     means the moment conditions cannot all hold at one value of the parameters.
#     It is a warning rather than an error because the estimates are the ones the
#     objective asked for: what is being reported on is the specification of the
#     system rather than the solve. Raised only where the two-step weight update
#     settled, since an unsettled weight matrix leaves J no reference distribution
#     to be judged against and the exhausted budget has been reported already.
#
#   deli_nested_solver_error
#     A rootSolve solve was asked for while another one was running. This one is
#     an error rather than a warning, because the solve cannot be attempted at
#     all rather than attempted and reported on. See `run_multiroot()`.

# Score magnitude below which a returned point is treated as solved, where the
# solver has not reported a convergence failure of its own. Both the
# iteration-budget check in solve_equations and the root-quality test in
# unsolved_equation compare against this floor, so it lives in one place. The
# floor is not an excuse on every path: see solve_equations() for why a solver
# that moved and then reported its own convergence test failing is judged on that
# report rather than on the size of the score it left behind.
score_floor <- 1e-4

# Largest share of an equation's per-observation contributions that may fail to
# cancel before the point is judged not to be a root. Applies only to equations
# whose contributions do not all carry one sign. See unsolved_equation().
noncancel_ceiling <- 0.9

# Largest value an equation whose contributions all carry one sign may take,
# measured against the scale of the problem, before the point is judged not to be
# a root. See unsolved_equation() for what the scale is and how the ceiling was
# chosen.
one_sided_ceiling <- 1e-4

# Largest Newton step a fit may still take at the point it returns, measured
# against the magnitude of the estimates themselves, before its equations are
# judged unsolved. Across the test suite, every fit that unsolved_equation()
# leaves quiet and whose bread can be solved takes a step either below 5e-4 or
# above 0.1, and the ceiling sits in that gap. The lower group is round-off and
# solver tolerance at points that are roots. The upper group opens with a
# logistic regression on separated data, whose slope runs away while its score
# falls away with it, so the contributions say nothing about it and this reading
# is the only one that sees it at all.
#
# The gap describes the fits this reading judges, not fits that stopped short
# generally. A fit unsolved_equation() catches first is never measured against
# the ceiling and its step bears no relation to it: the tobit fit of
# test-estimate-solvers.R is reported by its contributions at a step of 0.7. See
# relative_newton_step().
newton_step_ceiling <- 1e-3

# P-value below which the J-statistic of an over-identified GMM fit is reported as
# a rejection of its moment conditions. The chi-squared reference is asymptotic,
# so at a small sample and at a weight matrix that has not settled J runs
# over-sized, and the threshold has to clear that as well as catching a real
# failure. Measured over two hundred correctly specified fits at each size, the
# smallest P-value was 1.1e-4 at n = 15, 2.2e-5 at n = 30 and 2.9e-6 at n = 60,
# and 4.8e-6 at the identity weight matrix; the weakest signal from a badly
# specified fit was 8.8e-11. The threshold sits near the logarithmic middle of
# that gap, two orders of magnitude clear of the worst honest fit measured and two
# clear of the weakest real failure. See gmm_j_statistic().
j_rejection_p_value <- 1e-8

#' Find an equation that is not solved at a point
#'
#' The test cannot be the absolute size of the summed score, because that
#' depends on the scale of the estimating functions and on the number of
#' observations: a correctly solved stack whose per-observation contributions
#' run to 1e12 cannot drive its summed score below 1e-8 in double precision. Each
#' equation is therefore measured against a quantity carrying its own scale, and
#' which quantity that is depends on whether its contributions can cancel at all.
#' A score that is not finite is never a root.
#'
#' Where they can, the score is measured against the mass of the contributions
#' that produced it. The share of the mass that failed to cancel is near zero at
#' a genuine root, whatever the scale, and near one where the contributions all
#' point the same way.
#'
#' That share carries no information about an equation whose contributions all
#' carry one sign, because for such an equation it is exactly one at every point
#' except an exact zero, however close to the root the point is. Writing a
#' relation among the parameters is the way estimating equations reach that
#' state: see `vignette("custom-estimating-equations")` and the causal effect
#' rows of [ee_ipw()], [ee_aipw()] and [ee_gformula()], all of which repeat one
#' value across observations. Multiplying the same relation by an observation
#' weight leaves an equation that takes a different value at every observation
#' and still cannot cancel, so what matters is the one sign rather than the
#' repetition. What must vanish for such an equation is its mean contribution, so
#' that is what is measured, against the larger of the magnitude of the starting
#' values and the largest contribution anywhere in the stack. Taking the larger
#' of the two is what makes the reading free of scale: the first alone would grow
#' with the scale of the estimating functions, and the second alone is the value
#' itself for a stack that is nothing but one such equation, where the parameters
#' are the only other scale there is.
#'
#' The parameter term reads the starting values rather than the estimates because
#' the estimates are the thing under test. A solver that runs a parameter away
#' hands back a magnitude large enough to excuse whatever it left behind, so an
#' estimate-based scale is at its weakest in exactly the case this reading exists
#' to catch. The starting values are the caller's own statement of the scale the
#' parameters are expected on, and a fit that legitimately works at a magnitude
#' of 1e8 is started there, so the term keeps everything it was added for.
#'
#' Both readings are one-sided. Each is decisive when it is large, and when it is
#' small it has only failed to find evidence against the point: mixed-sign
#' contributions cancel well at points that are not roots at all, and a one-signed
#' equation whose value is negligible beside the rest of the stack is not judged
#' by this test however wrong the parameter it identifies. See
#' `relative_newton_step()` for what measures the rest.
#'
#' `one_sided_ceiling` was measured on the package's own fits. Across the test
#' suite, every article and vignette, and sweeps of the causal estimators and the
#' ratio stack over seeds, sample sizes and scalings, no correctly solved stack
#' carrying a one-signed equation read above 3.6e-8, while the fits that
#' genuinely stopped short read from 2.5e-4 to 1. The ceiling sits between them.
#'
#' @param ef A p-by-n matrix of per-observation contributions to the equations
#'   being solved, evaluated at the point under test.
#' @param theta The point under test.
#' @param init The starting values the solver was given, which supply the
#'   parameter term of the scale the one-sided reading measures against. Defaults
#'   to `theta`, for a caller judging a point that no solver produced.
#' @returns `NULL` when every equation is solved at that point. Otherwise a list
#'   naming the equation that misses its ceiling by the widest margin, in `row`,
#'   whether its summed score is finite, in `finite`, which of the two readings
#'   judged it, in `one_sided`, whether it has the same value at every
#'   observation, in `constant`, and the summed score it leaves, in `score`.
#' @noRd
unsolved_equation <- function(ef, theta, init = theta) {
  score <- rowSums(ef)
  if (!all(is.finite(score))) {
    bad <- which(!is.finite(score))[[1L]]
    return(list(
      row = bad,
      finite = FALSE,
      one_sided = FALSE,
      constant = FALSE,
      score = score[[bad]]
    ))
  }
  # Numerically negligible scores are always treated as solved, so equations
  # driven to zero (where neither reading means anything) never fail.
  if (max(abs(score)) <= score_floor) {
    return(NULL)
  }
  # An equation that cannot cancel is one with no negative contribution or none
  # positive.
  one_sided <- rowSums(ef < 0) == 0L | rowSums(ef > 0) == 0L
  # How far past its own ceiling each equation is, so the worst offender can be
  # named whichever of the two readings it failed.
  magnitude <- abs(ef)
  mass <- pmax(rowSums(magnitude), .Machine$double.eps)
  excess <- numeric(nrow(ef))
  excess[!one_sided] <- abs(score[!one_sided]) /
    mass[!one_sided] /
    noncancel_ceiling
  problem_scale <- max(1, abs(init[is.finite(init)]), magnitude)
  excess[one_sided] <- abs(score[one_sided]) /
    ncol(ef) /
    (one_sided_ceiling * problem_scale)
  worst <- which.max(excess)
  if (excess[[worst]] <= 1) {
    return(NULL)
  }
  list(
    row = worst,
    finite = TRUE,
    one_sided = one_sided[[worst]],
    constant = all(ef[worst, ] == ef[[worst, 1L]]),
    score = score[[worst]]
  )
}

#' Size of the Newton step towards a root of the estimating equations
#'
#' A just-identified problem has as many equations as parameters, so they must
#' vanish at a solution and the point a solver returns can be judged on how far
#' it still is from a root. Neither the size of the summed equations nor the
#' value of a minimizer's objective can say how far that is, because both carry
#' the scale of the estimating functions. The non-cancellation share used by
#' `unsolved_equation()` cannot say it either, for a different reason: where the
#' per-observation contributions are mixed-sign, which is every stack built from
#' a design matrix, that share stays near zero however wrong the parameters are.
#'
#' What says it is the bread, the derivative of the mean equations with respect
#' to the parameters, computed for the sandwich anyway. Solving it against the
#' mean equations gives the Newton step from the returned point, in the units of
#' the parameters themselves. Measuring each element against the magnitude of its
#' own estimate, with a floor of one so that an estimate which is legitimately
#' near zero cannot make a round-off step look large, leaves a number free of the
#' scale of the equations and of the scale of the parameters. At a solution it is
#' at round-off; where the solver stopped short of a root it is the distance
#' still to travel.
#'
#' A bread that cannot be solved leaves the point unjudged here. That is the
#' excuse owed to a non-differentiable estimating equation, whose Jacobian is
#' singular everywhere: where the derivative does not exist, the distance to a
#' root cannot be measured with it. The same excuse is made when a
#' finite-difference bread collapses to round-off, which is what happens when
#' `dx` is far smaller than the resolution of the contributions themselves.
#'
#' The excuse is owed to the equation that is solved at the point, not to the one
#' that is not. See `flat_equation()`, which judges the equations a singular
#' bread would otherwise carry past this reading unexamined.
#'
#' @param bread The bread matrix at the returned point, already divided by the
#'   number of observations.
#' @param moments The mean of each estimating equation at the returned point.
#' @param theta The returned parameter vector.
#' @returns The largest element of the Newton step, each measured against the
#'   magnitude of its own estimate with a floor of one, or `NA_real_` when the
#'   bread cannot be solved.
#' @noRd
relative_newton_step <- function(bread, moments, theta) {
  step <- tryCatch(solve(bread, moments), error = function(e) NULL)
  if (is.null(step) || !all(is.finite(step))) {
    return(NA_real_)
  }
  max(abs(step) / pmax(abs(theta), 1))
}

#' The summed value of each estimating equation
#'
#' The estimating function returns a p-by-n matrix, or a dimensionless vector
#' where there is a single equation, and both shapes have to sum the same way.
#'
#' @param vals An estimating-function return.
#' @returns A numeric vector with one entry per equation.
#' @noRd
equation_scores <- function(vals) {
  if (is.null(dim(vals))) sum(vals) else rowSums(vals)
}

#' Find an equation the bread has gone flat under and which is not solved
#'
#' A row of the bread that is identically zero says that the equation owning it
#' does not move when any parameter does, at the point under test and to the
#' resolution the derivative was taken at. Such a row makes the whole bread
#' singular, so `relative_newton_step()` cannot be taken and the point reaches
#' that reading only to be excused by it.
#'
#' The excuse is owed to an equation that is solved at the point, and to one
#' left where the caller put it. It is not owed to an equation that the solver
#' moved into a worse state than the one it was handed: a solver working from a
#' bread with this row in it has no reading of the equation at all, so any
#' movement it made was a trade against some other equation, and the value it
#' traded to is one nobody chose. No other reading catches that.
#' `unsolved_equation()` does not, because the contributions of the equations
#' this arises for are mixed-sign and cancel well at points that are not roots.
#' That is how a median equation left to a Levenberg-Marquardt solver comes back
#' at a value between two order statistics well away from the sample median,
#' with the deviation equation it was traded against solved to round-off.
#'
#' The comparison against the starting values is the same excuse
#' `solve_equations()` makes when rootSolve reports a failed convergence test
#' without having moved, and it is owed to the same estimating equations. A
#' median equation cannot be searched for at all, so the documentation of
#' [ee_positive_mean_deviation()] tells the caller to start `theta[2]` at the
#' sample median and expect the fit to hold there. Such a fit sums to a value
#' that is not zero whenever the number of observations is odd or the data are
#' tied, and it is the answer the caller asked for. The score floor is added to
#' the comparison so that a score reproduced to round-off counts as unchanged.
#'
#' The floor on the score itself is what keeps a row quiet that is flat because
#' its equation genuinely does not depend on the parameters. Such an equation is
#' at a root wherever it vanishes, and only the pair, a flat row and a score that
#' is not zero, is evidence of anything.
#'
#' @param ef A p-by-n matrix of per-observation contributions to the equations
#'   being solved, evaluated at the point under test.
#' @param bread The bread matrix at that point, with one row per equation.
#' @param init_score The summed value of each of those equations at the starting
#'   values, or `NULL` where the starting values are not known. Every judgement
#'   this reading makes rests on the comparison, so it reports nothing without
#'   them.
#' @returns `NULL` when no equation is flat, unsolved, and worse than it was at
#'   the starting values. Otherwise a list in the shape `unsolved_equation()`
#'   returns, naming the flat equation with the largest summed score in `row`,
#'   and carrying `flat = TRUE`.
#' @noRd
flat_equation <- function(ef, bread, init_score) {
  if (is.null(init_score)) {
    return(NULL)
  }
  zero_row <- vapply(
    seq_len(nrow(bread)),
    function(i) isTRUE(all(bread[i, ] == 0)),
    logical(1)
  )
  score <- equation_scores(ef)
  found <- which(
    zero_row &
      is.finite(score) &
      abs(score) > score_floor &
      abs(score) > abs(init_score) + score_floor
  )
  if (length(found) == 0L) {
    return(NULL)
  }
  worst <- found[[which.max(abs(score[found]))]]
  list(
    row = worst,
    finite = TRUE,
    one_sided = FALSE,
    constant = FALSE,
    flat = TRUE,
    score = score[[worst]],
    step = NA_real_
  )
}

#' Judge whether a returned point solves a just-identified system
#'
#' The three readings answer different questions and none subsumes another.
#' `unsolved_equation()` names an equation that cannot be at a root from its
#' contributions alone, which is the only thing that sees a stack whose Jacobian
#' does not exist. `relative_newton_step()` measures the distance still to
#' travel, which is the only thing that sees a stack built from a design matrix,
#' whose mixed-sign contributions cancel well however wrong the parameters are.
#' `flat_equation()` covers what falls between them: an equation whose
#' contributions cancel well at a point that is not its root, and whose own row
#' of the bread is what makes the Newton step unavailable.
#'
#' `flat_equation()` is asked before the Newton step because it explains why no
#' step can be taken, and asking it after would mean asking it only where the
#' step came back `NA` for that very reason.
#'
#' @param ef A p-by-n matrix of per-observation contributions to the equations
#'   being solved, evaluated at the point under test.
#' @param theta The point under test.
#' @param bread The bread matrix at that point, already divided by the number of
#'   observations.
#' @param init_score The summed value of each of those equations at the starting
#'   values, passed on to `flat_equation()`, which is the only reading that uses
#'   it. Defaults to `NULL`, which leaves that reading unmade.
#' @param init The starting values the solver was given, passed on to
#'   `unsolved_equation()`. Defaults to `theta`.
#' @returns `NULL` when the point solves the equations. Otherwise the list
#'   `unsolved_equation()` or `flat_equation()` returns, or, where neither
#'   reading found anything and the Newton step did, a list whose `row` is `NA`
#'   and whose `step` is that step. Every list carries `flat`, which is `TRUE`
#'   only for the second of the three.
#' @noRd
unsolved_point <- function(ef, theta, bread, init_score = NULL, init = theta) {
  found <- unsolved_equation(ef, theta, init)
  if (!is.null(found)) {
    found$flat <- FALSE
    found$step <- NA_real_
    return(found)
  }
  flat <- flat_equation(ef, bread, init_score)
  if (!is.null(flat)) {
    return(flat)
  }
  step <- relative_newton_step(bread, rowSums(ef) / ncol(ef), theta)
  if (!isTRUE(step > newton_step_ceiling)) {
    return(NULL)
  }
  list(
    row = NA_integer_,
    finite = TRUE,
    one_sided = FALSE,
    constant = FALSE,
    flat = FALSE,
    score = NA_real_,
    step = step
  )
}

#' Report that a solver returned a point that does not solve the equations
#'
#' The remedy differs by solver, and no message may send a user to rootSolve,
#' which is the solver that returns a spurious root silently.
#'
#' The diagnosis is reported alongside the solver's own account of itself only
#' where it says something the solver cannot. A flat equation is such a case: the
#' solver reports a convergence test that was met, since the equation it left
#' unsolved is invisible to the Jacobian it was working with, so the message has
#' to name the equation itself. The other readings say what the solver has
#' already said in its own terms, so they add no bullet.
#'
#' @param solved The list `solve_equations()` returned.
#' @param unsolved The list `unsolved_point()` returned, or `NULL` where the
#'   solver reported the failure itself and its own account of itself is what is
#'   reported.
#' @returns Invisible `NULL`, called for the warning.
#' @noRd
warn_unsolved <- function(solved, unsolved = NULL) {
  detail <- character(0)
  if (isTRUE(unsolved$flat)) {
    flat_row <- unsolved$row
    flat_score <- signif(unsolved$score, 3)
    detail <- c(
      "i" = "Estimating equation {flat_row} does not move when any parameter
             does, so no Newton step can measure the distance to its root, and
             it sums to {.val {flat_score}} rather than to zero."
    )
  }
  if (identical(solved$solver, "rootSolve")) {
    cli::cli_warn(
      c(
        "!" = "rootSolve did not converge to a root of the estimating
               equations.",
        "i" = "The estimating functions are not solved at the returned values
               (achieved precision {.val {signif(solved$precision, 3)}}).",
        detail,
        "i" = "Results may be unreliable. Consider using the {.val lm} solver or
               different starting values."
      ),
      class = "deli_solver_not_converged"
    )
  } else if (identical(solved$solver, "lm")) {
    cli::cli_warn(
      c(
        "!" = "minpack.lm stopped without solving the estimating equations
               (code {solved$code}).",
        "i" = "{solved$message}",
        detail,
        "i" = "Results may be unreliable. Consider increasing {.arg maxiter},
               using different starting values, or the {.val nleqslv} solver."
      ),
      class = "deli_solver_not_converged"
    )
  } else if (identical(solved$solver, "nleqslv")) {
    cli::cli_warn(
      c(
        "!" = "nleqslv stopped without solving the estimating equations
               (code {solved$code}).",
        "i" = "{solved$message}",
        detail,
        "i" = "Results may be unreliable. Consider using the {.val lm} solver or
               different starting values."
      ),
      class = "deli_solver_not_converged"
    )
  } else {
    cli::cli_warn(
      c(
        "!" = "The solver returned values that do not solve the estimating
               equations.",
        detail,
        "i" = "Results may be unreliable. Consider different starting values or
               one of the built-in solvers."
      ),
      class = "deli_solver_not_converged"
    )
  }
  invisible(NULL)
}

#' Report a solver's own convergence failure, worded by what the bread says
#'
#' `solve_equations()` takes `rootSolve::multiroot()` at its word when it reports
#' its convergence test failing after it moved, but that report describes the
#' search rather than the problem, and one thing that produces it is a problem
#' with no unique root to find. That diagnosis needs both halves of the evidence.
#' A rank-deficient bread says the mean estimating equations do not change along
#' at least one direction in the parameter space, so every point along it solves
#' them exactly as well as the returned one. A returned point that the readings
#' of `unsolved_point()` accept says the search is not what failed, since the
#' equations are satisfied where it stopped. Together they say the parameters are
#' not identified, and reporting a search that stopped short would then be true
#' of the search and beside the point about the problem.
#'
#' Both halves are needed because a bread can lose rank without the design being
#' the reason. Parameters that have run away leave a saturated estimating
#' function whose derivative underflows to zero, which is a whole column of
#' zeros; the point such a fit returns is not a root, and the solver's own
#' account of stopping short is the accurate one. So every failure but the
#' identified pair keeps that account.
#'
#' Neither wording can double up with the readings of the returned point.
#' `estimate_m_estimator()` reads the point once and reports it through exactly
#' one of the two branches, so one solve reports itself once.
#'
#' @param solved The list `solve_equations()` returned.
#' @param bread The bread matrix at the returned point, for the block of
#'   equations the solver worked on.
#' @param unsolved What `unsolved_point()` made of the returned point, `NULL`
#'   where it solves the equations.
#' @returns Invisible `NULL`, called for the warning.
#' @noRd
warn_solver_failure <- function(solved, bread, unsolved) {
  if (!is.null(unsolved) || !rank_deficient(bread)) {
    return(warn_unsolved(solved))
  }
  cli::cli_warn(
    c(
      "!" = "The estimating equations are rank deficient at the returned values,
             so the parameters are not identified.",
      "i" = "At least one direction in the parameter space leaves the mean
             estimating equations unchanged, so they have no unique root and
             every point along that direction solves them as well as the values
             returned here.",
      "i" = "Results may be unreliable. Check the equations for a redundant
             parameter and the design for linearly dependent columns."
    ),
    class = "deli_solver_not_converged"
  )
  invisible(NULL)
}

#' Whether a bread matrix leaves a direction the equations cannot see
#'
#' A square bread whose rank falls short of its dimension has at least one
#' direction in the parameter space along which the mean estimating equations do
#' not change, which is what a design with linearly dependent columns produces.
#'
#' A bread carrying values that are not finite is not read this way, both because
#' `qr()` cannot factor one and because a derivative that is not finite is a
#' different failure with a different remedy.
#'
#' @param bread The bread matrix at the returned point.
#' @returns `TRUE` where the matrix is square, finite, and rank deficient.
#' @noRd
rank_deficient <- function(bread) {
  square <- is.matrix(bread) && nrow(bread) == ncol(bread) && ncol(bread) > 0L
  if (!square || !all(is.finite(bread))) {
    return(FALSE)
  }
  qr(bread)$rank < ncol(bread)
}

#' Evaluate an expression with anything printed to stdout discarded
#'
#' Returns the value of `expr`, unlike [utils::capture.output()], which returns
#' the captured output and so forces the value to be assigned from inside the
#' call.
#'
#' @noRd
without_output <- function(expr) {
  null_con <- file(nullfile(), open = "wt")
  sink(null_con)
  on.exit(
    {
      sink()
      close(null_con)
    },
    add = TRUE
  )
  force(expr)
}

# ---- a solver's own account of itself ----------------------------------------
# A solver reports on itself through warnings as well as through its return
# value, and `solve_equations()` reads those reports and words deli's own. A
# report that reached the caller as well would describe one solve twice, in two
# vocabularies, so each is muffled where it is raised.
#
# What may not be muffled is a warning the estimating function raised. One fit
# evaluates that function at the starting values, at the point the solver
# returned, once or twice per parameter while the bread is built, and as many
# times as the solver likes in between. A warning from any of those is the
# caller's business, and the evaluations made inside the solver are no
# exception, so the muffling recognizes the solver's own reports by what they
# say and leaves everything else to propagate.
#
# The predicates below are what "recognizably the solver's own" means, one per
# solver that reports this way. They are vectorized over messages so that a
# recorded batch can be read in one call.

#' Whether `rootSolve::multiroot()` reported its convergence test failing
#'
#' The returned list carries no status flag, so the warning reading
#' "steady-state not reached" is the only account `multiroot()` gives of it. One
#' predicate decides both that the warning is the solver's own and so must not
#' reach the caller, and that the solve is to be reported as one whose
#' convergence test failed, so the two readings cannot come apart.
#'
#' @param messages A character vector of warning messages.
#' @returns A logical vector of the same length.
#' @noRd
multiroot_not_converged <- function(messages) {
  grepl("steady-state not reached", messages, fixed = TRUE)
}

#' Whether a message is `rootSolve::multiroot()` talking about itself
#'
#' Two reports exist. One is the convergence test above. The other comes from
#' the LINPACK routine that failed to factor the Jacobian, `dgefa` for a dense
#' matrix and `dgbfa` for a banded one, and its fixed text carries a run of
#' spaces before the reason. The pattern stops at the routine name rather than
#' matching through that run, so the number of spaces is not something this has
#' to know.
#'
#' @param messages A character vector of warning messages.
#' @returns A logical vector of the same length.
#' @noRd
multiroot_self_report <- function(messages) {
  multiroot_not_converged(messages) |
    grepl(
      "error during factorisation of matrix[[:space:]]*\\((dgefa|dgbfa)\\)",
      messages
    )
}

#' Whether a message is `minpack.lm::nls.lm()` talking about itself
#'
#' `nls.lm()` reports an exhausted budget twice: once as a bare warning out of
#' the MINPACK driver, carrying the info code, and once in the return value,
#' which is where `solve_equations()` reads it and words deli's own report. The
#' bare warning is written by the C code as the driver's name, the info code,
#' and the explanation, and the two drivers `nls.lm()` calls are `lmdif`, for a
#' Jacobian it approximates itself, and `lmder`, for one supplied to it.
#'
#' @param messages A character vector of warning messages.
#' @returns A logical vector of the same length.
#' @noRd
nls_lm_self_report <- function(messages) {
  grepl("^lm(dif|der): info = ", messages)
}

# ---- nested rootSolve solves -------------------------------------------------
# `rootSolve::multiroot()` cannot be called from inside itself. Its C code holds
# the environment it evaluates the function argument in in a single slot, so a
# second call started while the first is still running overwrites what the first
# is using, and the first then carries on over state belonging to the second.
# What that looks like depends on the two problems: the outer solve may fail out
# of the C code with a type error naming an environment, or it may return
# quietly with the inner fit's estimates in place of its own. The silent wrong
# answer is the worse of the two and neither is acceptable, so a rootSolve solve
# asked for while another one is running is refused instead of attempted.
#
# The marker below records that a `multiroot()` call is running, and it is set
# for the duration of that call and no longer. A fit evaluates its estimating
# function at the starting values, at the point the solver returned, and once or
# twice per parameter while the bread is built, all of them outside the solver;
# an inner rootSolve fit started from one of those is safe, and it is how a
# two-stage estimator whose outer fit runs on another solver is written.

solver_state <- new.env(parent = emptyenv())
solver_state$in_multiroot <- FALSE

#' Run `rootSolve::multiroot()` with the nested-solve marker set
#'
#' The marker is restored rather than cleared, and on the error path as well as
#' the ordinary one, so a solve that is refused or that fails leaves the marker
#' where it found it and the next fit solves normally. The restore is registered
#' before the marker is set, so that an interrupt arriving between the two
#' cannot leave the marker standing for the rest of the session; registering it
#' after would leave that window open.
#'
#' @param ... Arguments for [rootSolve::multiroot()].
#' @returns What `multiroot()` returns.
#' @noRd
run_multiroot <- function(...) {
  previous <- solver_state$in_multiroot
  on.exit(
    {
      solver_state$in_multiroot <- previous
    },
    add = TRUE
  )
  solver_state$in_multiroot <- TRUE
  rootSolve::multiroot(...)
}

#' Refuse a rootSolve solve asked for while another one is running
#'
#' @param call The frame to report the error against.
#' @returns Invisible `NULL` where no solve is running; otherwise it throws.
#' @noRd
check_not_nested_multiroot <- function(call) {
  if (!isTRUE(solver_state$in_multiroot)) {
    return(invisible(NULL))
  }
  cli::cli_abort(
    c(
      "!" = "A {.val rootSolve} solve cannot be nested inside another one.",
      "i" = "{.fn rootSolve::multiroot} evaluates the estimating function in a
             single environment it holds on to, so the inner solve overwrites
             what the outer one is still using and the outer one carries on
             over corrupted state.",
      "i" = "Give one of the two fits {.code solver = \"nleqslv\"}."
    ),
    class = "deli_nested_solver_error",
    call = call
  )
}

#' Evaluate a solver call with the solver's own warnings muffled and recorded
#'
#' @param expr The call to the solver.
#' @param self_report A predicate taking a character vector of warning messages
#'   and returning `TRUE` for each one that is the solver's own account of
#'   itself.
#' @returns A list holding the value of `expr` in `value` and the messages of
#'   the warnings that were muffled, in the order they were raised, in
#'   `reports`.
#' @noRd
with_muffled_self_reports <- function(expr, self_report) {
  reports <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      message <- conditionMessage(w)
      if (!self_report(message)) {
        # Not the solver's own, so it belongs to whoever wrote the estimating
        # function. Returning without muffling leaves it to the handlers
        # outside, which is how it reaches the caller.
        return()
      }
      reports <<- c(reports, message)
      # cnd_muffle() rather than invokeRestart("muffleWarning"), for the reason
      # given in R/conditions.R: a condition inheriting from "warning" that was
      # signalled without that restart would turn the restart call into an
      # error, while cnd_muffle() returns FALSE and leaves it to propagate.
      rlang::cnd_muffle(w)
    }
  )
  list(value = value, reports = reports)
}

#' Internal root-finding dispatcher
#'
#' A solver reports its own status in its own terms, and none of those terms is
#' a statement that the estimating equations are solved: every convergence test
#' on offer is relative to the residuals the solver started from, to the step it
#' last took, or to an absolute bound that the scale of the estimating functions
#' can put out of reach. So a report of failure is passed on as one, and
#' everything else is left for the caller to judge against the returned point
#' once the bread is in hand. See `unsolved_point()`.
#'
#' The one report taken at face value is `rootSolve::multiroot()` saying its
#' convergence test failed after it moved away from the starting values. A flat
#' estimating equation drives the summed score to nothing while the parameters
#' run away rather than settle, as a logistic regression on separated data does,
#' so at such a point neither reading of the returned values is evidence that
#' the equations are solved. The same report is not taken at face value when the
#' solver never moved, because a non-differentiable estimating equation, for
#' example the median in [ee_positive_mean_deviation()], has a singular Jacobian
#' everywhere and returns the starting values with that report and a summed score
#' that is not zero.
#'
#' Taking that report at face value is not the same as raising it here. What the
#' report means for the caller depends on the bread, which does not exist until
#' the solver has returned, so it is handed back in `not_converged` for the
#' caller to word. See `warn_solver_failure()`.
#'
#' Where a solver states its status in a warning rather than in its return
#' value, that warning is muffled here, because what it says is said again in
#' deli's own terms. Only what the solver says about itself is muffled: a
#' warning raised by `func` reaches the caller whether the evaluation that
#' raised it was made by this function or by the solver. See
#' `with_muffled_self_reports()`.
#'
#' @param func The summed estimating equations as a function of the parameters.
#' @param init The starting values.
#' @param method The name of the solver.
#' @param maxiter Iteration budget.
#' @param tolerance Solver tolerance.
#' @param call The frame the errors raised here report against, which is the
#'   [estimate()] method the user called rather than either of the internal
#'   functions between it and this one.
#' @returns A list holding the returned parameter vector in `par`, the name of
#'   the solver in `solver`, whether a warning has already been raised about the
#'   solve in `warned`, and whatever the solver reports that the warning needs:
#'   `code` and `message`, or `precision`. A `rootSolve` solve also carries
#'   `not_converged`, which is `TRUE` where the solver reported its own
#'   convergence test failing after moving and the caller has still to report it.
#' @noRd
solve_equations <- function(
  func,
  init,
  method,
  maxiter,
  tolerance,
  call = rlang::caller_env()
) {
  if (method == "rootSolve") {
    check_not_nested_multiroot(call)
    # rootSolve's Fortran code prints diagnostic messages to stdout, which stay
    # suppressed. Its warnings are read rather than suppressed wholesale: one of
    # them is the only report multiroot makes of its own convergence test
    # failing, since the returned list carries no status flag. If a future
    # rootSolve reworded it, the returned point still reaches the caller to be
    # judged rather than breaking. Warnings from `func` itself are the caller's
    # and pass through untouched; see `with_muffled_self_reports()`.
    run <- with_muffled_self_reports(
      without_output(
        run_multiroot(
          f = func,
          start = init,
          maxiter = maxiter,
          atol = tolerance
        )
      ),
      multiroot_self_report
    )
    result <- run$value
    not_converged <- any(multiroot_not_converged(run$reports))
    solved <- list(
      par = result$root,
      solver = "rootSolve",
      warned = FALSE,
      precision = result$estim.precis
    )
    # An exhausted iteration budget is a failure the solver can be believed
    # about: multiroot stops at a partially-solved point where the per-equation
    # contributions still cancel well, so the readings of the returned point
    # cannot see it, and the iteration count reaching maxiter is the reliable
    # signal. An exhausted budget is not a report that the convergence test
    # failed, so the score floor still excuses it.
    # isTRUE collapses a non-finite estim.precis (NaN or NA) to FALSE so the
    # comparison cannot leave an NA that would crash the if().
    budget_exhausted <- isTRUE(
      result$iter >= maxiter &&
        is.finite(result$estim.precis) &&
        result$estim.precis > score_floor
    )
    if (budget_exhausted) {
      cli::cli_warn(
        c(
          "!" = "rootSolve did not converge within {maxiter} iteration{?s}.",
          "i" = "The estimating functions are not solved at the returned values
                 (achieved precision {.val {signif(result$estim.precis, 3)}}).",
          "i" = "Results may be unreliable. Consider increasing {.arg maxiter}
                 or using the {.val lm} solver."
        ),
        class = "deli_solver_not_converged"
      )
      solved$warned <- TRUE
    } else {
      # The report is passed back rather than raised here, because what to say
      # about it depends on the bread and the bread does not exist yet. See
      # warn_solver_failure().
      solved$not_converged <- not_converged &&
        !isTRUE(all(result$root == init))
    }
    return(solved)
  }

  if (method == "lm") {
    rlang::check_installed("minpack.lm", reason = "for the lm solver.")
    # Python delicatessen solves with scipy.optimize.root(method = "lm"),
    # which drives MINPACK's lmdif. minpack.lm::nls.lm wraps the same MINPACK
    # code, minimising the sum of squares of the residual vector returned by
    # fn; at a root of the estimating equations that sum is zero. scipy passes
    # tol as xtol and maxiter as the function-evaluation budget (maxfev), so
    # ptol receives tolerance and maxfev receives maxiter. nls.lm caps its own
    # iteration counter at 1024, so the counter is set no higher to avoid a
    # spurious reset warning while maxfev carries the full budget.
    control <- minpack.lm::nls.lm.control(
      ptol = tolerance,
      maxfev = maxiter,
      maxiter = min(maxiter, 1024L)
    )
    # nls.lm reports an exhausted budget as a bare warning as well as in the
    # return value read below, so the bare one is muffled where it is raised.
    # Warnings from `func` itself pass through untouched.
    result <- with_muffled_self_reports(
      minpack.lm::nls.lm(par = init, fn = func, control = control),
      nls_lm_self_report
    )$value
    solved <- list(
      par = result$par,
      solver = "lm",
      warned = FALSE,
      code = result$info,
      message = result$message
    )
    # Codes 1 through 3 report a convergence test that was met, and each of the
    # three is relative to the residuals MINPACK started from, so a stack whose
    # contributions are large in magnitude meets one of them without moving.
    # Code 4 reports that the residual vector is orthogonal to every column of
    # the Jacobian to within gtol, which is where the algorithm can make no
    # further progress rather than where it has solved the system. None of the
    # four is evidence of a root, so all four leave the point to be judged.
    # Every other code is a failure of the algorithm itself and is reported as
    # one.
    if (!result$info %in% 1:4) {
      cli::cli_warn(
        c(
          "!" = "minpack.lm did not converge (code {result$info}).",
          "i" = "{result$message}",
          "i" = "Results may be unreliable. Consider increasing {.arg maxiter},
                 using different starting values, or the {.val nleqslv} solver."
        ),
        class = "deli_solver_not_converged"
      )
      solved$warned <- TRUE
    }
    return(solved)
  }

  if (method == "nleqslv") {
    rlang::check_installed("nleqslv", reason = "for the nleqslv solver.")
    result <- nleqslv::nleqslv(
      x = init,
      fn = func,
      control = list(maxit = maxiter, xtol = tolerance)
    )
    solved <- list(
      par = result$x,
      solver = "nleqslv",
      warned = FALSE,
      code = result$termcd,
      message = result$message
    )
    # termcd == 1 means nleqslv's function criterion was met. That criterion is
    # ftol, an absolute bound on the largest absolute function value with a
    # default of 1e-8, so a stack whose per-observation contributions are large
    # in magnitude cannot reach it however exactly the equations are solved, and
    # one whose contributions are small in magnitude meets it wherever it
    # happens to be. Such a fit stops on the x-value criterion and reports code
    # 2, or stalls and reports code 3. nleqslv documents both as leaving the
    # function values for the caller to judge rather than as failures, so codes
    # 1 to 3 all leave the point to be judged. Every other code is a failure of
    # the algorithm itself and is reported as one.
    #
    # Code 5 is a failure of the algorithm that says nothing against the point
    # it stopped at. nleqslv reports it where the Jacobian is too
    # ill-conditioned for it to take another step, and a stack whose blocks
    # differ by orders of magnitude reaches that state at the exact solution,
    # returning the values it started from to the last digit. What an
    # ill-conditioned Jacobian does put in doubt is the bread built there, and
    # so the variance, so that is what the report is about. It is still a
    # failure the solver has stated, so the returned point is not judged again
    # on top of it and the solve reports itself once.
    if (result$termcd == 5) {
      cli::cli_warn(
        c(
          "!" = "nleqslv reports an ill-conditioned Jacobian at the returned
                 values (code {result$termcd}).",
          "i" = "{result$message}",
          "i" = "The returned values may still solve the estimating equations;
                 it is the variance that an ill-conditioned Jacobian puts in
                 doubt. Consider rescaling the parameters or the estimating
                 functions."
        ),
        class = "deli_solver_not_converged"
      )
      solved$warned <- TRUE
    } else if (!result$termcd %in% 1:3) {
      cli::cli_warn(
        c(
          "!" = "nleqslv did not converge (code {result$termcd}).",
          "i" = "{result$message}",
          "i" = "Results may be unreliable. Consider using the {.val lm} solver
                 or different starting values."
        ),
        class = "deli_solver_not_converged"
      )
      solved$warned <- TRUE
    }
    return(solved)
  }

  cli::cli_abort("The solver {.val {method}} is not supported.", call = call)
}

# ---- GMMEstimator estimate method --------------------------------------------

#' @rdname estimate
#' @name estimate
method(estimate, GMMEstimator) <- function(
  object,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_dx(dx)
  # See the MEstimator method above and R/conditions.R for why the body sits in
  # a worker and what the scope does.
  without_repeated_warnings(estimate_gmm_estimator(
    object,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  ))
}

#' Minimize a `GMMEstimator` objective and assemble its sandwich variance
#'
#' The body of the [estimate()] method for [GMMEstimator()], separated from it so
#' that the method can wrap the whole of it in `without_repeated_warnings()`.
#'
#' @inheritParams estimate
#' @returns The `GMMEstimator` with its estimated properties populated.
#' @noRd
estimate_gmm_estimator <- function(
  object,
  solver,
  maxiter,
  tolerance,
  deriv_method,
  dx,
  allow_pinv
) {
  # Default solver for GMMEstimator is "BFGS" (optimization-based)
  if (is.null(solver)) {
    solver <- "BFGS"
  }
  stacked_equations <- object@stacked_equations
  init <- object@init
  subset <- object@subset
  finite_correction <- object@finite_correction

  # Evaluate stacked equations at init to determine dimensions, validating the
  # return before it reaches the objective and the sandwich components. The
  # validation also rejects an under-identified system. Its aborts name the
  # frame of the method that was called, as the MEstimator path does.
  vals_at_init <- eval_psi_at_init(
    stacked_equations,
    init,
    allow_over_identification = TRUE,
    error_call = rlang::caller_env()
  )
  if (is.null(dim(vals_at_init))) {
    # 1D case: single estimating equation
    n_eqs <- 1L
    n_obs <- length(vals_at_init)
  } else {
    n_eqs <- nrow(vals_at_init)
    n_obs <- ncol(vals_at_init)
  }

  # Determine if problem is over-identified
  over_identified <- length(init) < n_eqs

  # Initialize weight matrix as identity

  weight_matrix <- diag(n_eqs)

  # Current theta starts at init

  current_theta <- init

  # Build the GMM objective function
  # Minimizes: t(sum_ef / n) %*% Q %*% (sum_ef / n)
  build_gmm_objective <- function(weight_mat, current_full_theta) {
    function(theta) {
      # If subset, expand theta into full parameter vector
      if (!is.null(subset)) {
        full_theta <- current_full_theta
        full_theta[subset] <- theta
      } else {
        full_theta <- theta
      }

      # Evaluate estimating equations
      ef <- stacked_equations(full_theta)
      if (is.null(dim(ef))) {
        sum_ef <- sum(ef)
      } else {
        sum_ef <- rowSums(ef)
      }

      # Scale by n for numerical stability (as Python does)
      sum_ef <- sum_ef / n_obs

      # GMM objective: t(sum_ef) %*% Q %*% sum_ef
      as.numeric(crossprod(sum_ef, weight_mat %*% sum_ef))
    }
  }

  # Get initial values for optimization (possibly subset)
  if (!is.null(subset)) {
    inits <- current_theta[subset]
  } else {
    inits <- current_theta
  }

  # STEP 1.1: Initial minimization with identity weight matrix
  # The refusal minimize_gmm() raises is about the solver the caller asked for,
  # so it names the frame of the method that was called rather than either of
  # the internal functions in between, as the MEstimator path does.
  gmm_obj <- build_gmm_objective(weight_matrix, current_theta)
  solved <- minimize_gmm(
    gmm_obj,
    inits,
    solver,
    maxiter,
    tolerance,
    call = rlang::caller_env()
  )
  theta_solved <- solved$par
  minimizer_converged <- solved$converged

  # Expand theta if subset was used
  if (!is.null(subset)) {
    current_theta[subset] <- theta_solved
  } else {
    current_theta <- theta_solved
  }

  # STEP 1.2: Iterative procedure for over-identified problems
  #
  # Whether the update settled is what decides later on whether the J-statistic
  # has a reference distribution to be judged against. A budget of zero updates
  # leaves the identity weight matrix in place without ever reporting a failure,
  # so it counts as settled here and its J is judged too; the threshold is set
  # clear of what a correctly specified fit reports at the identity.
  weight_update_settled <- TRUE
  if (over_identified) {
    overid_maxiter <- object@overid_maxiter
    overid_tolerance <- object@overid_tolerance

    if (overid_maxiter >= 1L) {
      for (iter in seq_len(overid_maxiter)) {
        prev_theta <- current_theta

        # Update weight matrix: Q = inverse of meat matrix
        evald_q <- stacked_equations(current_theta)
        if (is.null(dim(evald_q))) {
          # Reshape a vector-return psi to 1-by-n so the meat cross-product is
          # 1-by-1, mirroring the MEstimator path.
          evald_q <- matrix(evald_q, nrow = 1)
        }
        meat_q <- compute_meat(evald_q) / n_obs
        if (allow_pinv) {
          weight_matrix <- tryCatch(
            solve(meat_q),
            error = function(e) {
              rlang::check_installed(
                "MASS",
                reason = "for pseudo-inverse when the meat matrix is singular."
              )
              MASS::ginv(meat_q)
            }
          )
        } else {
          weight_matrix <- solve(meat_q)
        }

        # Re-minimize with updated weight matrix
        if (!is.null(subset)) {
          inits <- current_theta[subset]
        } else {
          inits <- current_theta
        }

        gmm_obj <- build_gmm_objective(weight_matrix, current_theta)
        solved <- minimize_gmm(
          gmm_obj,
          inits,
          solver,
          maxiter,
          tolerance,
          call = rlang::caller_env()
        )
        theta_solved <- solved$par
        minimizer_converged <- solved$converged

        # Expand theta if subset was used
        if (!is.null(subset)) {
          current_theta[subset] <- theta_solved
        } else {
          current_theta <- theta_solved
        }

        # Check convergence
        error <- max(abs(current_theta - prev_theta))
        if (error <= overid_tolerance) {
          break
        }

        # Warn if max iterations reached without convergence
        if (iter == overid_maxiter && error > overid_tolerance) {
          weight_update_settled <- FALSE
          cli::cli_warn(c(
            "!" = "{.arg overid_maxiter} ({overid_maxiter}) has been exceeded
                   for the iterative GMM updating.",
            "i" = "Terminating the over-identification procedure at the last
                   estimated values."
          ))
        }
      }
    }
  }

  # STEP 2: Compute sandwich variance
  evald <- stacked_equations(current_theta)
  if (is.null(dim(evald))) {
    # Reshape a vector-return psi to 1-by-n so the meat cross-product is
    # 1-by-1, mirroring the MEstimator path.
    evald <- matrix(evald, nrow = 1)
  }

  # Bread matrix
  bread <- compute_bread(stacked_equations, current_theta, deriv_method, dx) /
    n_obs

  # A just-identified problem has as many moment conditions as parameters, so
  # the moments must vanish at a solution. minimize_gmm only inspects optim's
  # status code, and optim reports convergence on the flat tail of an objective
  # that never reaches zero, so nothing so far has looked at the moments.
  #
  # A minimizer that reported a failure of its own has already warned about it,
  # so its point is not judged again. An over-identified problem cannot drive
  # every moment to zero, and a subset fit holds the parameters outside the
  # subset at their initial values while still summing every equation into the
  # objective. None of those three can be judged this way, so none is.
  if (minimizer_converged && !over_identified && is.null(subset)) {
    unsolved <- unsolved_point(
      evald,
      current_theta,
      bread,
      init_score = equation_scores(vals_at_init),
      init = init
    )
    if (!is.null(unsolved)) {
      summed_moment <- signif(unsolved$score, 3)
      detail <- if (unsolved$flat) {
        # A flat row has a row number, so this branch comes ahead of the Newton
        # step, and its contributions say nothing, so it comes ahead of the
        # readings of them too.
        "Moment condition {unsolved$row} does not move when any parameter does,
         so no Newton step can measure the distance to its root, and it sums to
         {.val {summed_moment}} rather than to zero."
      } else if (is.na(unsolved$row)) {
        "A Newton step from the estimated values would move at least one of them
         by {.val {signif(unsolved$step, 3)}} times the larger of its own
         magnitude and one, so the minimizer stopped short of a root."
      } else if (!unsolved$finite) {
        "The summed value of moment condition {unsolved$row} is
         {.val {summed_moment}}, which is not finite."
      } else if (unsolved$constant) {
        # A constant row is one-signed too, so this branch has to come first. It
        # earns its own wording because repeating one value is how the vignettes
        # and the causal estimators write a relation among the parameters.
        "Moment condition {unsolved$row} has the same value at every
         observation and does not vanish, leaving a summed moment of
         {.val {summed_moment}}."
      } else if (unsolved$one_sided) {
        "The per-observation contributions of moment condition {unsolved$row}
         all carry one sign, so they cannot cancel, and their mean does not
         vanish, leaving a summed moment of {.val {summed_moment}}."
      } else {
        "The per-observation contributions of moment condition {unsolved$row}
         do not cancel, leaving a summed moment of {.val {summed_moment}}."
      }
      cli::cli_warn(
        c(
          "!" = "The moment conditions are not solved at the estimated values.",
          "i" = detail,
          "i" = "Results may be unreliable. Consider different starting values
                 or another {.arg solver}."
        ),
        class = "deli_solver_not_converged"
      )
    }
  }

  # An over-identified problem cannot be judged that way, because no value of the
  # parameters drives every moment to zero and a residual moment is expected
  # rather than diagnostic. Hansen's J-statistic is the reading it has instead:
  # the size of what is left over, measured against a reference distribution
  # rather than against the scale of the data. See GMMEstimator() for what it
  # means and where its reference distribution holds.
  #
  # A subset fit is left unjudged for the reason the check above leaves it
  # unjudged: the parameters outside the subset were held at their initial values
  # rather than estimated, and the degrees of freedom the reference distribution
  # takes assume every parameter was.
  j_statistic <- NULL
  if (over_identified && is.null(subset)) {
    j_statistic <- gmm_j_statistic(evald, weight_matrix, n_obs)
    if (weight_update_settled) {
      j_df <- n_eqs - length(current_theta)
      j_p_value <- stats::pchisq(j_statistic, j_df, lower.tail = FALSE)
      if (j_p_value < j_rejection_p_value) {
        cli::cli_warn(
          c(
            "!" = "The over-identifying moment conditions are rejected at the
                   estimated values.",
            "i" = "Hansen's J-statistic is {.val {signif(j_statistic, 4)}} on
                   {j_df} degree{?s} of freedom, for a P-value of
                   {.val {signif(j_p_value, 3)}}.",
            "i" = "A J-statistic this large usually means the moment conditions
                   cannot all hold at one value of the parameters, so it is the
                   conditions that need checking rather than the starting
                   values."
          ),
          class = "deli_gmm_moments_rejected"
        )
      }
    }
  }

  # Meat matrix
  meat_mat <- compute_meat(evald) / n_obs
  meat_mat <- finite_sample_correction(
    meat_mat,
    n_obs,
    length(current_theta),
    finite_correction
  )

  # Sandwich variance
  asymp_var <- build_sandwich(bread, meat_mat, allow_pinv)

  if (is.null(asymp_var)) {
    var_mat <- NULL
  } else {
    var_mat <- asymp_var / n_obs
  }

  # Apply parameter names, from `init` where it has any and otherwise from the
  # row names of the evaluation already made above. No second call to the
  # estimating functions is needed to read them. A refused fill names the frame
  # of the method that was called, as the MEstimator path does.
  param_names <- resolve_param_names(
    names(object@init),
    evald,
    length(current_theta),
    error_call = rlang::caller_env()
  )
  names(current_theta) <- param_names
  # Bread and meat may be non-square in the over-identified case
  # (n_eqs x n_params for bread, n_eqs x n_eqs for meat). Where they are, they
  # are indexed by the moment conditions rather than by the parameters, so they
  # are left unlabeled rather than carrying labels a reader could line up
  # against coef(). Row names on the estimating functions used to reach the meat
  # regardless, through tcrossprod(), and the weight matrix through the solve()
  # of it; the strips below are what keep the three matrices consistent.
  if (nrow(bread) == ncol(bread) && nrow(bread) == length(param_names)) {
    dimnames(bread) <- list(param_names, param_names)
  } else {
    dimnames(bread) <- NULL
  }
  if (
    nrow(meat_mat) == ncol(meat_mat) && nrow(meat_mat) == length(param_names)
  ) {
    dimnames(meat_mat) <- list(param_names, param_names)
  } else {
    dimnames(meat_mat) <- NULL
  }
  # The weight matrix is indexed by the moment conditions on every path, and is
  # the identity on a just-identified fit, so it is never parameter-labeled.
  dimnames(weight_matrix) <- NULL
  if (!is.null(asymp_var)) {
    dimnames(asymp_var) <- list(param_names, param_names)
  }
  if (!is.null(var_mat)) {
    dimnames(var_mat) <- list(param_names, param_names)
  }

  # Update the object
  object@theta <- current_theta
  object@n_obs <- as.integer(n_obs)
  object@bread <- bread
  object@meat <- meat_mat
  object@asymptotic_variance <- asymp_var
  object@variance <- var_mat
  object@weight_matrix <- weight_matrix
  object@j_statistic <- j_statistic

  object
}

#' Hansen's J-statistic at a GMM minimum
#'
#' n times the GMM objective at the estimated values, which is the quantity the
#' chi-squared reference distribution of an over-identified fit judges. See
#' [GMMEstimator()] for what it reports and where that distribution holds.
#'
#' @param evald The estimating functions at the estimated values, as a matrix of
#'   one row per moment condition.
#' @param weight_matrix The weight matrix the fit finished with.
#' @param n_obs The number of observations.
#'
#' @returns The statistic, as a single number.
#' @noRd
gmm_j_statistic <- function(evald, weight_matrix, n_obs) {
  # The same averaged moments the objective is built from, so the statistic is n
  # times the value the minimizer reported rather than a second reading of it.
  moments <- rowSums(evald) / n_obs
  n_obs * as.numeric(crossprod(moments, weight_matrix %*% moments))
}

#' Internal minimization dispatcher for GMMEstimator
#'
#' @param call The frame to report the refusal of an unusable `solver` against.
#'   That refusal is about what the caller asked for, and this dispatcher is two
#'   frames below the method the caller reached, so the frame travels with the
#'   arguments.
#'
#' @returns A list holding the minimizing parameter vector in `par` and, in
#'   `converged`, whether the minimizer reported success. A minimizer that
#'   reported a failure has already warned about it, so the caller leaves the
#'   moments at the point it returned unjudged.
#' @noRd
minimize_gmm <- function(func, init, method, maxiter, tolerance, call) {
  if (is.character(method)) {
    # Use stats::optim for minimization
    result <- stats::optim(
      par = init,
      fn = func,
      method = method,
      control = list(maxit = maxiter, reltol = tolerance)
    )
    if (result$convergence != 0) {
      cli::cli_warn(
        c(
          "!" = "Minimization did not converge (code {result$convergence}).",
          "i" = "{result$message}"
        ),
        class = "deli_solver_not_converged"
      )
    }
    return(list(par = result$par, converged = result$convergence == 0))
  }

  if (is.function(method)) {
    result <- method(stacked_equations = func, init = init)
    check_solver_return(result, length(init))
    # A custom minimizer reports nothing, so its point is judged like one a
    # built-in minimizer declared solved.
    return(list(par = result, converged = TRUE))
  }

  cli::cli_abort("{.arg solver} must be a string or function.", call = call)
}

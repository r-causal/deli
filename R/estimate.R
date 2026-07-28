#' Estimate parameters and sandwich variance
#'
#' Solves the estimating equations for the parameter vector `theta` and
#' computes the empirical sandwich variance estimator.
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
#'   fits are exceptions: an over-identified one cannot drive every moment to
#'   zero, and a `subset` one holds some parameters fixed while still summing
#'   every equation into the objective, so neither is judged either way and
#'   neither warns. For those, inspect `rowSums()` of the estimating functions
#'   at the returned values.
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
#'   `deriv_method = "exact"` (default 1e-9).
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
#' m <- MEstimator(stacked_equations = psi, init = c(0)) |>
#'   estimate()
#' m@theta
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
  # error instead of an opaque failure inside the solver.
  check_psi_at_init(stacked_equations(init), init)

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
    solved <- solve_equations(summed_ee, inits, solver, maxiter, tolerance)
  } else if (is.function(solver)) {
    par <- solver(stacked_equations = summed_ee, init = inits)
    check_solver_return(par, length(inits))
    # A custom solver reports nothing, so its point is judged like one a
    # built-in solver declared solved.
    solved <- list(par = par, solver = "custom", warned = FALSE)
  } else {
    cli::cli_abort("{.arg solver} must be a string or function.")
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

  # Judge the returned point, unless the solver reported a failure of its own
  # and has already warned about it. This waits for the bread because the
  # Jacobian is what measures the distance still to travel, and a solver whose
  # own convergence test is relative to its starting residuals can report
  # success without moving at all. A subset fit solves only the subset
  # equations for the subset parameters, so only that block is judged.
  if (!solved$warned) {
    judged <- subset %||% seq_len(nrow(evald))
    unsolved <- unsolved_point(
      evald[judged, , drop = FALSE],
      full_theta[judged],
      bread[judged, judged, drop = FALSE]
    )
    if (!is.null(unsolved)) {
      warn_unsolved(solved)
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

  # Apply parameter names
  param_names <- names(object@init) %||% paste0("theta_", seq_along(full_theta))
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

# ---- solver convergence conditions -------------------------------------------
# The warnings raised when a solver does not solve the estimating equations
# carry a condition class, so a caller or a test can match on the class rather
# than on the prose, as the exact-mode aborts in `R/autodiff.R` do.
#
#   deli_solver_not_converged
#     The solver stopped without solving the estimating equations, either
#     because it reported a failure of its own or because the returned point
#     does not solve them. Every solver branch raises this one class. A solver
#     that reported a failure warns from solve_equations() and its point is not
#     judged again, so one M-estimation solve raises the class at most once. A
#     just-identified GMM solve judges its moments only where the minimizer
#     reported success, so it too raises it at most once; an over-identified
#     solve calls the minimizer repeatedly and can raise it once per pass.

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
# above 20, and the ceiling sits in that gap. The gap describes the fits this
# reading judges, not fits that stopped short generally: a fit
# unsolved_equation() catches first is never measured against the ceiling and
# can sit well below it, the smallest in the suite being the tobit fit of
# test-estimate-solvers.R at 0.7. See relative_newton_step().
newton_step_ceiling <- 1

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
#' that is what is measured, against the larger of the magnitude of the estimates
#' and the largest contribution anywhere in the stack. Taking the larger of the
#' two is what makes the reading free of scale: the first alone would grow with
#' the scale of the estimating functions, and the second alone is the value
#' itself for a stack that is nothing but one such equation, where the estimates
#' are the only scale available.
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
#' @returns `NULL` when every equation is solved at that point. Otherwise a list
#'   naming the equation that misses its ceiling by the widest margin, in `row`,
#'   whether its summed score is finite, in `finite`, which of the two readings
#'   judged it, in `one_sided`, whether it has the same value at every
#'   observation, in `constant`, and the summed score it leaves, in `score`.
#' @keywords internal
#' @noRd
unsolved_equation <- function(ef, theta) {
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
  problem_scale <- max(1, abs(theta[is.finite(theta)]), magnitude)
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
#' A bread that cannot be solved leaves the point unjudged. That is the excuse
#' owed to a non-differentiable estimating equation, whose Jacobian is singular
#' everywhere: where the derivative does not exist, the distance to a root cannot
#' be measured with it. The same excuse is made when a finite-difference bread
#' collapses to round-off, which is what happens when `dx` is far smaller than
#' the resolution of the contributions themselves.
#'
#' @param bread The bread matrix at the returned point, already divided by the
#'   number of observations.
#' @param moments The mean of each estimating equation at the returned point.
#' @param theta The returned parameter vector.
#' @returns The largest element of the Newton step, each measured against the
#'   magnitude of its own estimate with a floor of one, or `NA_real_` when the
#'   bread cannot be solved.
#' @keywords internal
#' @noRd
relative_newton_step <- function(bread, moments, theta) {
  step <- tryCatch(solve(bread, moments), error = function(e) NULL)
  if (is.null(step) || !all(is.finite(step))) {
    return(NA_real_)
  }
  max(abs(step) / pmax(abs(theta), 1))
}

#' Judge whether a returned point solves a just-identified system
#'
#' The two readings answer different questions and neither subsumes the other.
#' `unsolved_equation()` names an equation that cannot be at a root from its
#' contributions alone, which is the only thing that sees a stack whose Jacobian
#' does not exist. `relative_newton_step()` measures the distance still to
#' travel, which is the only thing that sees a stack built from a design matrix,
#' whose mixed-sign contributions cancel well however wrong the parameters are.
#'
#' @param ef A p-by-n matrix of per-observation contributions to the equations
#'   being solved, evaluated at the point under test.
#' @param theta The point under test.
#' @param bread The bread matrix at that point, already divided by the number of
#'   observations.
#' @returns `NULL` when the point solves the equations. Otherwise the list
#'   `unsolved_equation()` returns, or, where that reading found nothing and the
#'   Newton step did, a list whose `row` is `NA` and whose `step` is that step.
#' @keywords internal
#' @noRd
unsolved_point <- function(ef, theta, bread) {
  found <- unsolved_equation(ef, theta)
  if (!is.null(found)) {
    found$step <- NA_real_
    return(found)
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
    score = NA_real_,
    step = step
  )
}

#' Report that a solver returned a point that does not solve the equations
#'
#' The remedy differs by solver, and no message may send a user to rootSolve,
#' which is the solver that returns a spurious root silently.
#'
#' @param solved The list `solve_equations()` returned.
#' @returns Invisible `NULL`, called for the warning.
#' @keywords internal
#' @noRd
warn_unsolved <- function(solved) {
  if (identical(solved$solver, "rootSolve")) {
    cli::cli_warn(
      c(
        "!" = "rootSolve did not converge to a root of the estimating
               equations.",
        "i" = "The estimating functions are not solved at the returned values
               (achieved precision {.val {signif(solved$precision, 3)}}).",
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
        "i" = "Results may be unreliable. Consider different starting values or
               one of the built-in solvers."
      ),
      class = "deli_solver_not_converged"
    )
  }
  invisible(NULL)
}

#' Evaluate an expression with anything printed to stdout discarded
#'
#' Returns the value of `expr`, unlike [utils::capture.output()], which returns
#' the captured output and so forces the value to be assigned from inside the
#' call.
#'
#' @keywords internal
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
#' @param func The summed estimating equations as a function of the parameters.
#' @param init The starting values.
#' @param method The name of the solver.
#' @param maxiter Iteration budget.
#' @param tolerance Solver tolerance.
#' @returns A list holding the returned parameter vector in `par`, the name of
#'   the solver in `solver`, whether a warning has already been raised about the
#'   solve in `warned`, and whatever the solver reports that the warning needs:
#'   `code` and `message`, or `precision`.
#' @keywords internal
#' @noRd
solve_equations <- function(func, init, method, maxiter, tolerance) {
  if (method == "rootSolve") {
    # rootSolve's Fortran code prints diagnostic messages to stdout, and the
    # call may emit warnings of its own as well as any raised inside func.
    # Everything stays suppressed, as before, but the messages are kept, because
    # one of them is the only report multiroot makes of its own convergence
    # test failing: the returned list carries no status flag. If a future
    # rootSolve reworded it, the returned point still reaches the caller to be
    # judged rather than breaking.
    solver_messages <- character()
    result <- withCallingHandlers(
      without_output(
        rootSolve::multiroot(
          f = func,
          start = init,
          maxiter = maxiter,
          atol = tolerance
        )
      ),
      warning = function(w) {
        solver_messages <<- c(solver_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    not_converged <- any(grepl(
      "steady-state not reached",
      solver_messages,
      fixed = TRUE
    ))
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
    } else if (not_converged && !isTRUE(all(result$root == init))) {
      warn_unsolved(solved)
      solved$warned <- TRUE
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
    result <- minpack.lm::nls.lm(par = init, fn = func, control = control)
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
    if (!result$termcd %in% 1:3) {
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

  cli::cli_abort("The solver {.val {method}} is not supported.")
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
  # Default solver for GMMEstimator is "BFGS" (optimization-based)
  if (is.null(solver)) {
    solver <- "BFGS"
  }
  stacked_equations <- object@stacked_equations
  init <- object@init
  subset <- object@subset
  finite_correction <- object@finite_correction

  # Evaluate stacked equations at init to determine dimensions, validating the
  # return before it reaches the objective and the sandwich components.
  vals_at_init <- stacked_equations(init)
  check_psi_at_init(vals_at_init, init, allow_over_identification = TRUE)
  if (is.null(dim(vals_at_init))) {
    # 1D case: single estimating equation
    n_eqs <- 1L
    n_obs <- length(vals_at_init)
  } else {
    n_eqs <- nrow(vals_at_init)
    n_obs <- ncol(vals_at_init)
  }

  # Check for under-identification
  n_init <- length(init)
  if (n_init > n_eqs) {
    cli::cli_abort(
      "The number of initial values ({n_init}) must be less than or equal to
       the number of estimating equations ({n_eqs})."
    )
  }

  # Determine if problem is over-identified
  over_identified <- n_init < n_eqs

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
  gmm_obj <- build_gmm_objective(weight_matrix, current_theta)
  solved <- minimize_gmm(gmm_obj, inits, solver, maxiter, tolerance)
  theta_solved <- solved$par
  minimizer_converged <- solved$converged

  # Expand theta if subset was used
  if (!is.null(subset)) {
    current_theta[subset] <- theta_solved
  } else {
    current_theta <- theta_solved
  }

  # STEP 1.2: Iterative procedure for over-identified problems
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
        solved <- minimize_gmm(gmm_obj, inits, solver, maxiter, tolerance)
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
    unsolved <- unsolved_point(evald, current_theta, bread)
    if (!is.null(unsolved)) {
      summed_moment <- signif(unsolved$score, 3)
      detail <- if (is.na(unsolved$row)) {
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

  # Apply parameter names
  param_names <- names(object@init) %||%
    paste0("theta_", seq_along(current_theta))
  names(current_theta) <- param_names
  # Bread and meat may be non-square in the over-identified case
  # (n_eqs x n_params for bread, n_eqs x n_eqs for meat)
  if (nrow(bread) == ncol(bread) && nrow(bread) == length(param_names)) {
    dimnames(bread) <- list(param_names, param_names)
  }
  if (
    nrow(meat_mat) == ncol(meat_mat) && nrow(meat_mat) == length(param_names)
  ) {
    dimnames(meat_mat) <- list(param_names, param_names)
  }
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

  object
}

#' Internal minimization dispatcher for GMMEstimator
#'
#' @returns A list holding the minimizing parameter vector in `par` and, in
#'   `converged`, whether the minimizer reported success. A minimizer that
#'   reported a failure has already warned about it, so the caller leaves the
#'   moments at the point it returned unjudged.
#' @keywords internal
#' @noRd
minimize_gmm <- function(func, init, method, maxiter, tolerance) {
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

  cli::cli_abort("{.arg solver} must be a string or function.")
}

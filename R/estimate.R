#' Estimate parameters and sandwich variance
#'
#' Solves the estimating equations for the parameter vector `theta` and
#' computes the empirical sandwich variance estimator.
#'
#' @param object An `MEstimator` object.
#' @param solver Character string specifying the solver algorithm, or a custom
#'   function. When `NULL` (default), uses `"rootSolve"` for `MEstimator`
#'   ([rootSolve::multiroot()]) and `"BFGS"` for `GMMEstimator`
#'   ([stats::optim()]). Other options for `MEstimator`: `"lm"`, the
#'   Levenberg-Marquardt algorithm ([minpack.lm::nls.lm()]), which mirrors the
#'   default solver of Python `delicatessen`
#'   (`scipy.optimize.root(method = "lm")`); and `"nleqslv"` (uses
#'   [nleqslv::nleqslv()]). A custom function must accept
#'   `stacked_equations` and `init` arguments and return the solved theta
#'   vector.
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
#' @param ... Additional arguments (currently unused).
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
#' m <- MEstimator(stacked_equations = psi, init = c(0))
#' m <- estimate(m)
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

  # Per-observation contributions of the subset equations at a parameter
  # vector, used to judge whether the returned root actually solves the
  # equations (see solve_equations()).
  matrix_ee <- function(theta) {
    if (!is.null(subset)) {
      full_theta <- init
      full_theta[subset] <- theta
    } else {
      full_theta <- theta
    }
    ef <- stacked_equations(full_theta)
    if (is.null(dim(ef))) {
      ef <- matrix(ef, nrow = 1)
    }
    if (!is.null(subset)) {
      ef <- ef[subset, , drop = FALSE]
    }
    ef
  }

  # Get initial values (possibly subset)
  if (!is.null(subset)) {
    inits <- init[subset]
  } else {
    inits <- init
  }

  # Solve
  if (is.character(solver)) {
    theta_solved <- solve_equations(
      summed_ee,
      inits,
      solver,
      maxiter,
      tolerance,
      ee_matrix = matrix_ee
    )
  } else if (is.function(solver)) {
    theta_solved <- solver(stacked_equations = summed_ee, init = inits)
    check_solver_return(theta_solved, length(inits))
  } else {
    cli::cli_abort("{.arg solver} must be a string or function.")
  }

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

#' Internal root-finding dispatcher
#' @keywords internal
#' @noRd
solve_equations <- function(
  func,
  init,
  method,
  maxiter,
  tolerance,
  ee_matrix = NULL
) {
  if (method == "rootSolve") {
    # rootSolve's Fortran code prints diagnostic messages to stdout
    # and may emit warnings such as a singular matrix. Suppress both.
    suppressWarnings(invisible(capture.output({
      result <- rootSolve::multiroot(
        f = func,
        start = init,
        maxiter = maxiter,
        atol = tolerance
      )
    })))
    # rootSolve::multiroot reports success even when it does not reach a root,
    # and the suppression above discards its own "steady-state not reached"
    # warning. Two distinct failure modes are surfaced here, matching the lm and
    # nleqslv branches. First, an exhausted iteration budget: multiroot stops at
    # a partially-solved point where the per-equation contributions still cancel
    # well, so the cancellation heuristic below cannot see it; the iteration
    # count reaching maxiter is the reliable signal.
    score_floor <- 1e-4
    if (result$iter >= maxiter && result$estim.precis > score_floor) {
      cli::cli_warn(c(
        "!" = "rootSolve did not converge within {maxiter} iteration{?s}.",
        "i" = "The estimating functions are not solved at the returned values
               (achieved precision {.val {signif(result$estim.precis, 3)}}).",
        "i" = "Results may be unreliable. Consider increasing {.arg maxiter} or
               using the {.val lm} solver."
      ))
    } else if (!is.null(ee_matrix)) {
      # Second, a singular Jacobian at the starting values: multiroot returns
      # those values unchanged with a non-zero score. Inspect the score at the
      # returned root and warn when the equations are not solved there.
      warn_if_not_root(result, ee_matrix)
    }
    return(result$root)
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
    # info values 1 through 4 indicate a successful convergence test
    if (!result$info %in% 1:4) {
      cli::cli_warn(c(
        "!" = "minpack.lm did not converge (code {result$info}).",
        "i" = "{result$message}",
        "i" = "Results may be unreliable. Consider using the default
               {.val rootSolve} solver instead."
      ))
    }
    return(result$par)
  }

  if (method == "nleqslv") {
    rlang::check_installed("nleqslv", reason = "for the nleqslv solver.")
    result <- nleqslv::nleqslv(
      x = init,
      fn = func,
      control = list(maxit = maxiter, xtol = tolerance)
    )
    # termcd == 1 means converged; anything else is a problem
    if (result$termcd != 1) {
      cli::cli_warn(c(
        "!" = "nleqslv did not converge (code {result$termcd}).",
        "i" = "{result$message}",
        "i" = "Results may be unreliable. Consider using the default
               {.val rootSolve} solver instead."
      ))
    }
    return(result$x)
  }

  cli::cli_abort("The solver {.val {method}} is not supported.")
}

#' Warn when rootSolve returns a point that does not solve the equations
#'
#' `rootSolve::multiroot()` reports success even when it makes no progress from
#' the starting values, which happens when the estimating-function Jacobian is
#' singular there. The returned score is then far from zero and the reported
#' `root` is spurious. This helper inspects the returned point and warns.
#'
#' The test cannot be the absolute size of the summed score, because it depends
#' on the scale and number of observations, and some non-differentiable
#' estimating equations (for example the median in
#' [ee_positive_mean_deviation()]) legitimately settle at a point where the
#' summed score is not zero. Instead it measures, for each equation, how much
#' of the per-observation contributions cancel: at a genuine root they cancel
#' almost completely, whereas a solver that never left a degenerate point
#' leaves contributions that all point the same way. The warning fires only
#' when the summed score is not numerically negligible and the largest
#' non-cancellation fraction is close to one.
#'
#' @param result The list returned by [rootSolve::multiroot()].
#' @param ee_matrix A function returning the per-observation matrix of the
#'   solved equations at a parameter vector.
#' @keywords internal
#' @noRd
warn_if_not_root <- function(result, ee_matrix) {
  score <- result$f.root
  # Numerically negligible scores are always treated as solved, so equations
  # driven to zero (where the cancellation fraction is meaningless) never warn.
  score_floor <- 1e-4
  if (max(abs(score)) <= score_floor) {
    return(invisible(NULL))
  }
  ef <- ee_matrix(result$root)
  mass <- rowSums(abs(ef))
  # Fraction of each equation's contributions that failed to cancel. Near one
  # means no cancellation occurred, the signature of a non-root.
  noncancel <- abs(rowSums(ef)) / pmax(mass, .Machine$double.eps)
  if (max(noncancel) > 0.9) {
    cli::cli_warn(c(
      "!" = "rootSolve did not converge to a root of the estimating equations.",
      "i" = "The estimating functions are not solved at the returned values
             (achieved precision {.val {signif(result$estim.precis, 3)}}).",
      "i" = "Results may be unreliable. Consider using the {.val lm} solver or
             different starting values."
    ))
  }
  invisible(NULL)
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
  theta_solved <- minimize_gmm(gmm_obj, inits, solver, maxiter, tolerance)

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
        theta_solved <- minimize_gmm(gmm_obj, inits, solver, maxiter, tolerance)

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
      cli::cli_warn(c(
        "!" = "Minimization did not converge (code {result$convergence}).",
        "i" = "{result$message}"
      ))
    }
    return(result$par)
  }

  if (is.function(method)) {
    result <- method(stacked_equations = func, init = init)
    check_solver_return(result, length(init))
    return(result)
  }

  cli::cli_abort("{.arg solver} must be a string or function.")
}

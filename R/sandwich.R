#' Compute the bread matrix
#'
#' Computes the bread matrix for the empirical sandwich variance estimator.
#' The bread is the negative Jacobian of the summed estimating equations,
#' scaled by `1/n`.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions.
#' @param theta Numeric vector of parameter estimates.
#' @param deriv_method Character string for the derivative method. One of
#'   `"capprox"` (central), `"fapprox"` (forward), or `"bapprox"` (backward).
#' @param dx Numeric step size (default `1e-9`). The step is absolute, floored
#'   at the floating-point resolution of each estimate; see
#'   [approx_differentiation()].
#'
#' @returns A p-by-p bread matrix.
#'
#' @keywords internal
compute_bread <- function(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-9
) {
  deriv_method <- check_deriv_method(deriv_method)

  # Sum the estimating equations across observations
  summed_ee <- function(input_theta) {
    ef <- stacked_equations(input_theta)
    # Handle PrimalTangent returns (from autodiff-compatible EE functions)
    if (is_pt(ef)) {
      return(sum(ef))
    }
    # A tangent-carrying matrix (the usual p-by-n return of a built-in EE under
    # the exact pass): sum across observations within each equation, carrying
    # the tangents in parallel. Checked before the list branch because a
    # PrimalTangentArray is structurally a list whose `[[` yields scalar pairs.
    #
    # A dimensionless tangent array (a psi built with `c()`, or any vector
    # return) sums every element into one value, whereas the list branch below
    # sums within each element. The two agree only when there is a single
    # equation, which is the only case either shape can reach: `estimate()`
    # treats a vector return as one 1-by-n equation, so a p > 1 fit whose psi
    # returns a vector fails on the resulting 1-by-p Jacobian regardless.
    if (is_pt_array(ef)) {
      if (is.null(dim(ef$primal))) {
        return(primal_tangent(sum(ef$primal), sum(ef$tangent)))
      }
      return(primal_tangent_array(rowSums(ef$primal), rowSums(ef$tangent)))
    }
    if (is.list(ef) && any(vapply(ef, is_pt, logical(1)))) {
      # A plain list holding one value per equation, from an `lapply()` or
      # `sapply()` inside the psi. `c()` returns a PrimalTangentArray and so takes
      # the branch above instead.
      #
      # Every element is reduced with `sum()`, which each tangent surface has a
      # method for, so an equation that arrives as a scalar pair, as a tangent
      # array, or as a plain numeric constant all sum across observations here.
      # A scalar pair anywhere in the list is what says the list is per-equation:
      # reading the first element alone decided that for the whole list, so a psi
      # whose theta-free equation or whose `c()`-built equation came first was
      # refused the reduction this branch carries out.
      return(lapply(ef, sum))
    }
    # Under the exact pass, any other list shape is one this step cannot sum.
    # `base::rbind()` on a single tangent-carrying value, which is what a psi
    # written in the global environment reaches, builds a 1-by-2 list matrix
    # holding one slot per cell. Neither cell holds a pair, which is why it misses
    # the per-equation branch above; its `dim` is why it would go on to reach
    # `rowSums()` below, which fails with the opaque `'x' must be numeric or
    # complex`. Route it to a diagnostic instead.
    #
    # The two ways a list can arrive here are different failures. A list whose
    # elements still carry tangents, such as an `lapply()` returning one tangent
    # array per equation, has lost nothing and only its container shape is
    # unsupported, so it gets its own condition rather than being told the
    # derivatives are gone.
    if (deriv_method == "exact" && is.list(ef)) {
      if (any(vapply(ef, is_tangent_container, logical(1)))) {
        pt_unsupported_shape_abort()
      }
      pt_tangent_lost_abort()
    }
    if (is.null(dim(ef))) {
      return(sum(ef))
    }
    rowSums(ef)
  }

  if (deriv_method == "exact") {
    bread_matrix <- auto_differentiation(theta, summed_ee)
  } else {
    bread_matrix <- approx_differentiation(
      func = summed_ee,
      theta = theta,
      method = deriv_method,
      dx = dx
    )
  }

  if (anyNA(bread_matrix)) {
    cli::cli_warn(c(
      "!" = "The bread matrix contains NA values, so it cannot be inverted.",
      "i" = "The variance will not be calculated."
    ))
  }

  -1 * bread_matrix
}

#' Compute the meat matrix
#'
#' Computes the meat matrix as the cross-product of the estimating equation
#' evaluations: \eqn{EE \times EE^T}.
#'
#' @param evaluations A p-by-n matrix of estimating equation evaluations, where
#'   p is the number of parameters and n is the number of observations.
#'
#' @returns A p-by-p meat matrix.
#'
#' @keywords internal
compute_meat <- function(evaluations) {
  tcrossprod(evaluations)
}

#' Build the sandwich variance estimator
#'
#' Combines bread and meat matrices into the sandwich:
#' \eqn{B^{-1} M (B^{-1})^T}.
#'
#' @param bread A p-by-p bread matrix.
#' @param meat A p-by-p meat matrix.
#' @param allow_pinv Logical. If `TRUE` (default), uses the pseudo-inverse
#'   when the bread matrix is singular.
#'
#' @returns A p-by-p sandwich covariance matrix, or `NULL` if the bread
#'   contains `NA` values.
#'
#' @keywords internal
build_sandwich <- function(bread, meat, allow_pinv = TRUE) {
  # If bread contains NA, return NULL

  if (anyNA(bread)) {
    return(NULL)
  }

  # Invert the bread matrix
  if (allow_pinv) {
    # Use pseudo-inverse for robustness
    bread_inv <- tryCatch(
      solve(bread),
      error = function(e) {
        rlang::check_installed(
          "MASS",
          reason = "for pseudo-inverse when the bread matrix is singular."
        )
        MASS::ginv(bread)
      }
    )
  } else {
    bread_inv <- solve(bread)
  }

  # Sandwich: B^{-1} M (B^{-1})^T
  bread_inv %*% meat %*% t(bread_inv)
}

#' Apply finite-sample correction to the meat matrix
#'
#' Applies the HC1 correction: \eqn{meat \times n / (n - p)}.
#'
#' @param meat A p-by-p meat matrix.
#' @param n Integer number of observations.
#' @param p Integer number of parameters.
#' @param adjustment Character string or `NULL`. Currently only `"HC1"` is
#'   supported.
#'
#' @returns The corrected meat matrix.
#'
#' @keywords internal
finite_sample_correction <- function(meat, n, p, adjustment = NULL) {
  if (is.null(adjustment)) {
    return(meat)
  }

  adj_upper <- toupper(adjustment)
  if (adj_upper == "HC1") {
    if (n <= p) {
      cli::cli_abort(
        "The number of observations ({n}) is not greater than the number of
         parameters ({p}), so the HC1 correction cannot be applied."
      )
    }
    meat * n / (n - p)
  } else {
    cli::cli_abort(
      c(
        "The finite-sample correction {.val {adjustment}} is not available.",
        "i" = "Supported options: {.val NULL}, {.val HC1}."
      )
    )
  }
}

#' Compute the empirical sandwich variance estimator
#'
#' Computes the empirical sandwich variance estimator directly from a set of
#' estimating equations and a vector of parameter estimates. Unlike
#' [MEstimator()], this function does not solve for the parameters; it assumes
#' that `theta` is already the root of the estimating equations and only
#' assembles the sandwich covariance at that point.
#'
#' The sandwich is built from a bread matrix and a meat matrix. The bread is the
#' negative Jacobian of the summed estimating equations, and the meat is the
#' cross-product of the equation evaluations. Each is scaled by `1/n` internally,
#' so the returned matrix is on the asymptotic scale (see `Value`). The bread
#' Jacobian is obtained either by finite differences or by forward-mode
#' automatic differentiation.
#'
#' @param stacked_equations A function that takes a numeric vector `theta` and
#'   returns a p-by-n matrix of estimating equation contributions, where p is
#'   the number of parameters and n is the number of observations.
#' @param theta Numeric vector of parameter estimates. This function assumes
#'   `theta` is the root of `stacked_equations`; it does not solve for it.
#' @param deriv_method Character string selecting the method used to build the
#'   bread Jacobian. One of `"capprox"` (central difference), `"fapprox"`
#'   (forward difference), `"bapprox"` (backward difference), or `"exact"`
#'   (forward-mode automatic differentiation). Default `"capprox"`.
#'
#'   This default differs from Python delicatessen, whose `compute_sandwich`
#'   defaults to `"approx"`, a forward difference computed through SciPy's
#'   `approx_fprime`. deli does not replicate the SciPy `"approx"` path; its
#'   `"fapprox"` is the hand-implemented forward difference. Code ported across
#'   the two libraries should set `deriv_method` explicitly rather than rely on
#'   the default.
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. A small value is recommended, since large steps
#'   can give poor approximations. Default `1e-9`. Must be a single positive
#'   finite number, which is checked whichever `deriv_method` is in force. The
#'   step is absolute and is floored at the floating-point resolution of each
#'   estimate, so a large parameter magnitude cannot silently reduce it to
#'   nothing; see [approx_differentiation()].
#' @param allow_pinv Logical. When `TRUE` (default), the Moore-Penrose
#'   pseudo-inverse is used if the bread matrix is singular; when `FALSE`, a
#'   singular bread matrix raises an error.
#' @param finite_correction Character string or `NULL`. Finite-sample correction
#'   applied to the meat matrix. `NULL` (default) applies no correction; `"HC1"`
#'   rescales the meat by \eqn{n / (n - p)}, where p is the number of parameters.
#'
#' @returns A p-by-p covariance matrix on the asymptotic scale. The bread and
#'   meat are each divided by n internally, so the returned matrix is the
#'   variance that corresponds to the standard deviation. Dividing it by the
#'   number of observations gives the standard-error-scale variance, whose
#'   square-rooted diagonal is the vector of standard errors.
#'
#' @seealso [MEstimator()] and [GMMEstimator()], which solve for `theta` and
#'   report this variance internally, and [delta_method()] for the variance of a
#'   transformation of the parameters.
#'
#' @examples
#' # A generic data set for estimating a mean and variance
#' y <- c(1, 2, 4, 1, 2, 3, 1, 5, 2)
#'
#' # The mean and variance are the roots of ee_mean_variance, so they can be
#' # computed directly rather than solved for
#' theta <- c(mean(y), stats::var(y) * (length(y) - 1) / length(y))
#'
#' # Wrap the built-in estimating equation as a function of theta alone
#' psi <- function(theta) ee_mean_variance(theta, y = y)
#'
#' # compute_sandwich() returns the asymptotic-scale variance, so dividing by
#' # n puts it on the standard-error scale
#' sandwich <- compute_sandwich(psi, theta = theta) / length(y)
#' sandwich
#'
#' # The diagonal square roots are the standard errors
#' sqrt(diag(sandwich))
#'
#' @references
#' Boos DD, & Stefanski LA. (2013). M-estimation (estimating equations). In
#' Essential Statistical Inference (pp. 297-337). Springer, New York, NY.
#'
#' @export
compute_sandwich <- function(
  stacked_equations,
  theta,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  finite_correction = NULL
) {
  check_dx(dx)
  # This evaluates the estimating function once for itself and once or twice per
  # parameter for the bread, so an estimating function that warns raises the same
  # warning several times for one call. See R/conditions.R. The body is short
  # enough to wrap in place, unlike the estimate() methods, which put theirs in a
  # worker.
  without_repeated_warnings({
    # Evaluate estimating equations at theta-hat
    evald <- stacked_equations(theta)
    if (is.null(dim(evald))) {
      n_obs <- length(evald)
      n_params <- 1
      # Reshape to 1-by-n matrix so tcrossprod works correctly in compute_meat
      evald <- matrix(evald, nrow = 1)
    } else {
      n_obs <- ncol(evald)
      n_params <- nrow(evald)
    }

    # Step 1: Bread matrix
    bread <- compute_bread(stacked_equations, theta, deriv_method, dx)
    bread <- bread / n_obs

    # Step 2: Meat matrix
    meat <- compute_meat(evald)
    meat <- meat / n_obs
    meat <- finite_sample_correction(meat, n_obs, n_params, finite_correction)

    # Step 3: Build sandwich
    build_sandwich(bread, meat, allow_pinv)
  })
}

#' Confidence bands for parameter vectors
#'
#' Computes simultaneous confidence bands that provide coverage for the
#' entire parameter vector, adjusting for multiple comparisons. The formula is:
#' \deqn{\hat{\theta} \pm \hat{c}_{\alpha/2} \times \widehat{SE}(\hat{\theta})}
#' where \eqn{\hat{c}} is the adjusted critical value.
#'
#' @param object A fitted `MEstimator` object, or a numeric vector of
#'   parameter estimates.
#' @param alpha Numeric significance level, between 0 and 1. Default `0.05`.
#' @param method Character string. `"supt"` (default) for supremum-t or
#'   `"bonferroni"` for Bonferroni correction.
#' @param n_draws Integer number of MVN draws for the sup-t method.
#'   Default `1e5`. The critical value is a Monte Carlo estimate of a fixed
#'   sup-t quantile, and at `1e5` draws the band half-width varies by roughly
#'   0.2% across independent draws, which is negligible relative to the band
#'   width. Python's estimator method defaults to `1e6`; both defaults target
#'   the same quantity and agree to within Monte Carlo error, so `deli` keeps
#'   the smaller default for faster computation. Set `n_draws = 1e6` to match
#'   Python's estimator default. Note that seeded band values are not
#'   reproducible across the two languages because the MVN samplers differ.
#' @param seed Integer seed for reproducibility. Default `NULL`.
#' @param subset Integer vector of parameter indices to compute bands for.
#'   Default `NULL` (all parameters).
#' @param covariance Numeric covariance matrix (only when `object` is numeric).
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored.
#'
#' @returns A p-by-2 matrix with columns `"lower"` and `"upper"`. For a fitted
#'   estimator the rows are named for the parameters, as in
#'   [`confint()`][deli-generics]. For a numeric `object` they take their names
#'   from it, when it has any.
#'
#' @examplesIf requireNamespace("MASS", quietly = TRUE)
#' # Two independent samples, each contributing one mean parameter
#' set.seed(42)
#' n <- 200
#' y1 <- rnorm(n, 2)
#' y2 <- rnorm(n, 3)
#'
#' psi <- function(theta) {
#'   rbind(y1 - theta[1], y2 - theta[2])
#' }
#'
#' m <- m_estimate(stacked_equations = psi, init = c(0, 0))
#'
#' # Simultaneous bands cover both means at once, so they are wider than the
#' # pointwise intervals from confint(m)
#' confidence_bands(m, method = "supt", seed = 1)
#'
#' @export
# n_draws defaults to 1e5, not the 1e6 Python uses on its estimator method. The
# sup-t critical value converges to the same quantile either way, and at 1e5 the
# residual Monte Carlo error is negligible relative to the band width while
# running about an order of magnitude faster. Set n_draws = 1e6 for numerical
# agreement with Python under matched settings.
confidence_bands <- new_generic(
  "confidence_bands",
  "object",
  function(
    object,
    alpha = 0.05,
    method = "supt",
    n_draws = 100000L,
    seed = NULL,
    subset = NULL,
    covariance = NULL,
    ...
  ) {
    S7::S7_dispatch()
  }
)

method(confidence_bands, deli_estimator) <- function(
  object,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL,
  subset = NULL,
  covariance = NULL,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)

  if (!is.null(subset)) {
    theta <- object@theta[subset]
    covariance <- object@variance[subset, subset, drop = FALSE]
  } else {
    theta <- object@theta
    covariance <- object@variance
  }

  compute_confidence_bands(
    theta,
    covariance,
    alpha = alpha,
    method = method,
    n_draws = n_draws,
    seed = seed
  )
}

method(confidence_bands, class_numeric) <- function(
  object,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL,
  subset = NULL,
  covariance = NULL,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  if (!is.null(subset)) {
    object <- object[subset]
    if (!is.null(covariance)) {
      covariance <- as.matrix(covariance)[subset, subset, drop = FALSE]
    }
  }

  compute_confidence_bands(
    object,
    covariance,
    alpha = alpha,
    method = method,
    n_draws = n_draws,
    seed = seed
  )
}

#' Compute confidence bands from theta and covariance
#'
#' @param theta Numeric parameter vector.
#' @param covariance Numeric covariance matrix.
#' @param alpha Significance level. Default `0.05`.
#' @param method `"supt"` or `"bonferroni"`. Default `"supt"`.
#' @param n_draws Number of MVN draws for sup-t. Default `1e5`. See
#'   [confidence_bands()] for why this differs from Python's `1e6` estimator
#'   default and how to match it.
#' @param seed RNG seed. Default `NULL`.
#'
#' @returns A p-by-2 matrix with columns `"lower"` and `"upper"`. Rows take
#'   their names from `theta`, when it has any.
#'
#' @examplesIf requireNamespace("MASS", quietly = TRUE)
#' fit <- m_estimate(mpg ~ wt + hp, data = mtcars, .ee = ee_regression,
#'                   model = "linear")
#'
#' # Bands from the estimates and covariance alone, without the fitted object
#' compute_confidence_bands(coef(fit), covariance = vcov(fit),
#'                          method = "supt", seed = 1)
#'
#' @export
compute_confidence_bands <- function(
  theta,
  covariance,
  alpha = 0.05,
  method = "supt",
  n_draws = 100000L,
  seed = NULL
) {
  # as.numeric() strips the names, so they are kept aside for the row labels of
  # the returned bands, which are the same labels confint() puts on its rows.
  param_names <- names(theta)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)
  k <- length(theta)

  # Validate dimensions and alpha before computing, mirroring Python's
  # confidence_bands. Without these, a mismatched theta recycles over the
  # covariance and an out-of-range alpha returns infinite or NaN bands.
  n_rows <- nrow(covariance)
  n_cols <- ncol(covariance)
  if (n_rows != n_cols) {
    cli::cli_abort(
      "{.arg covariance} must be a square matrix, but it has {n_rows} row{?s}
       and {n_cols} column{?s}."
    )
  }
  if (k != n_rows) {
    cli::cli_abort(
      "{.arg theta} and {.arg covariance} must share a dimension, but
       {.arg theta} has length {k} and {.arg covariance} has dimension
       {n_rows}."
    )
  }
  check_alpha(alpha)

  se <- sqrt(diag(covariance))

  method <- tolower(method)
  if (method %in% c("supt", "sup-t")) {
    if (any(se <= 0)) {
      cli::cli_abort(
        "At least one parameter has a standard error of zero or less. The sup-t method cannot be applied."
      )
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }
    rlang::check_installed(
      "MASS",
      reason = "to draw from the multivariate normal for sup-t confidence bands."
    )
    # Draw from multivariate normal
    draws <- MASS::mvrnorm(n = n_draws, mu = rep(0, k), Sigma = covariance)
    if (k == 1) {
      draws <- matrix(draws, ncol = 1)
    }
    # Scale by SE and take absolute value
    scaled <- abs(sweep(draws, 2, se, `/`))
    # Max across parameters for each draw
    ts <- apply(scaled, 1, max)
    # Critical value is (1-alpha) quantile
    cv <- quantile(ts, probs = 1 - alpha, names = FALSE)
  } else if (method == "bonferroni") {
    cv <- qnorm(1 - alpha / (2 * k))
  } else {
    cli::cli_abort(
      "Method {.val {method}} is not supported. Use {.val supt} or {.val bonferroni}."
    )
  }

  cb <- cbind(lower = theta - cv * se, upper = theta + cv * se)
  rownames(cb) <- param_names
  cb
}

#' Delta method for variance of transformed parameters
#'
#' Computes the variance-covariance matrix for a transformation of parameters
#' using the Delta Method:
#' \deqn{Var[g(\theta)] \approx G \Sigma G^T}
#' where \eqn{G} is the Jacobian of \eqn{g} and \eqn{\Sigma} is the
#' covariance matrix of \eqn{\theta}.
#'
#' @param object A fitted `MEstimator` object, or a numeric vector of
#'   parameter estimates.
#' @param transform Function that takes `theta` and returns a numeric vector.
#' @param covariance Numeric covariance matrix (only used when `object` is a
#'   numeric vector).
#' @param deriv_method Character string for the derivative method used to build
#'   the Jacobian of `transform`. One of `"capprox"` (central difference),
#'   `"fapprox"` (forward difference), `"bapprox"` (backward difference), or
#'   `"exact"` (forward-mode automatic differentiation). Default `"capprox"`.
#' @param dx Numeric step size for the finite-difference methods; ignored when
#'   `deriv_method = "exact"`. Default `1e-9`. Must be a single positive finite
#'   number, which is checked whichever `deriv_method` is in force. The step is
#'   absolute and is floored at the floating-point resolution of each estimate,
#'   so a large parameter magnitude cannot silently reduce it to nothing; see
#'   [approx_differentiation()].
#' @param ... Not used. Must be empty, so a name that is not one of the
#'   documented arguments is an error rather than silently ignored. Exact names
#'   matter here because `deriv_method` selects how the Jacobian is built, and a
#'   dropped misspelling would leave the default in place and return a different
#'   variance with nothing to signal the substitution.
#'
#' @returns A covariance matrix for `transform(theta)`.
#'
#' @examples
#' fit <- m_estimate(vs ~ mpg, data = mtcars, .ee = ee_regression,
#'                   model = "logistic")
#'
#' # Variance of the odds ratio for mpg, exponentiating the log-odds coefficient
#' delta_method(fit, transform = function(theta) exp(theta[2]))
#'
#' # The same variance from the estimates and covariance alone
#' delta_method(coef(fit), transform = function(theta) exp(theta[2]),
#'              covariance = vcov(fit))
#'
#' @export
delta_method <- new_generic(
  "delta_method",
  "object",
  function(
    object,
    transform,
    covariance = NULL,
    deriv_method = "capprox",
    dx = 1e-9,
    ...
  ) {
    S7::S7_dispatch()
  }
)

method(delta_method, deli_estimator) <- function(
  object,
  transform,
  covariance = NULL,
  deriv_method = "capprox",
  dx = 1e-9,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  check_estimated(object)
  # One call differentiates `transform` and so evaluates it several times, once
  # for the shape check and once or twice per parameter for the Jacobian. A
  # transform that warns would otherwise repeat itself for one call. See
  # R/conditions.R. The scope wraps the shared worker rather than this body for
  # the same reason it does in the estimate() methods.
  without_repeated_warnings(
    delta_method_impl(
      object@theta,
      transform,
      object@variance,
      deriv_method,
      dx
    )
  )
}

method(delta_method, class_numeric) <- function(
  object,
  transform,
  covariance = NULL,
  deriv_method = "capprox",
  dx = 1e-9,
  ...
) {
  rlang::check_dots_empty(call = rlang::caller_env())
  # See the method above for why the worker call sits inside the scope.
  without_repeated_warnings(
    delta_method_impl(object, transform, covariance, deriv_method, dx)
  )
}

#' @noRd
delta_method_impl <- function(
  theta,
  g,
  covariance,
  deriv_method = "capprox",
  dx = 1e-9
) {
  # Validate covariance presence before any coercion. as.matrix(NULL) otherwise
  # fails with an opaque "'data' must be of a vector type" error before the
  # reachable shape checks below.
  if (is.null(covariance)) {
    cli::cli_abort(c(
      "{.arg covariance} is required but was not supplied.",
      "i" = "Pass the covariance matrix of the parameter estimates when calling
             {.fn delta_method} on a numeric vector."
    ))
  }
  deriv_method <- check_deriv_method(deriv_method)
  check_dx(dx)
  theta <- as.numeric(theta)
  covariance <- as.matrix(covariance)

  # Validate shapes before computing, mirroring the checks in Python's
  # delta_method. Without these, a bad shape either returns a wrongly shaped
  # matrix or fails with an opaque non-conformable error.
  theta_star <- g(theta)
  star_dim <- dim(theta_star)
  # A plain numeric vector has no dim. R's matrix algebra yields column vectors
  # (n x 1), which are one-dimensional in intent, so those are allowed; a
  # genuine matrix (more than one column) or higher-dimensional array is not.
  if (!is.null(star_dim) && (length(star_dim) > 2 || star_dim[2] > 1)) {
    cli::cli_abort(
      "Output from {.arg transform} must be a one-dimensional vector."
    )
  }
  n_rows <- nrow(covariance)
  n_cols <- ncol(covariance)
  if (n_rows != n_cols) {
    cli::cli_abort(
      "{.arg covariance} must be a square matrix, but it has {n_rows} row{?s}
       and {n_cols} column{?s}."
    )
  }
  v <- length(theta)
  if (n_rows != v) {
    cli::cli_abort(
      "{.arg covariance} and the parameter vector must share a dimension, but
       the parameter vector has length {v} and {.arg covariance} has dimension
       {n_rows}."
    )
  }

  # Delta method: G * Sigma * G^T
  g_prime <- transform_jacobian(theta, g, deriv_method, dx)
  g_prime %*% covariance %*% t(g_prime)
}

#' The Jacobian of a transform at a parameter vector
#'
#' The derivative half of the delta method, separated from the product that
#' follows it because not every caller wants the whole covariance matrix. A
#' caller predicting one quantity per row of a design needs the variance of each
#' row on its own, which is the diagonal of that product, and forming the whole
#' matrix to take its diagonal costs the square of the number of rows in memory.
#' Such a caller takes the Jacobian from here and forms
#' `rowSums((G %*% Sigma) * G)` instead.
#'
#' `deriv_method` is expected already normalized by [check_deriv_method()].
#'
#' @param theta The parameter vector to differentiate at.
#' @param g The transform.
#' @param deriv_method The normalized derivative method.
#' @param dx The finite-difference step.
#' @returns A matrix with one row per element of `g(theta)` and one column per
#'   parameter.
#' @noRd
transform_jacobian <- function(theta, g, deriv_method, dx) {
  # Exact autodiff must receive the raw transform: as.numeric() would strip the
  # tangents and yield a zero Jacobian. extract_tangent_column() already handles
  # every shape the raw transform can return. The finite-difference methods need
  # plain numeric output, so they see the wrapped one.
  if (deriv_method == "exact") {
    return(auto_differentiation(theta, g))
  }
  approx_differentiation(
    func = function(t) as.numeric(g(t)),
    theta = theta,
    method = deriv_method,
    dx = dx
  )
}

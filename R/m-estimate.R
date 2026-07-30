#' One-step M-estimation
#'
#' Creates an `MEstimator` and estimates it in one call, analogous to how
#' [stats::lm()] creates and fits a model in a single step. Supports both a
#' formula interface (for regression-family estimating equations) and a
#' function interface (for custom estimating equations).
#'
#' @details
#' Both interfaces place `...` ahead of `init`, `subset`, `tolerance`, and the
#' rest of the settings, so each of those must be named in full: R does not
#' partially match a supplied name against an argument that follows `...`.
#' [estimate()] and the inference generics take `...` last, so they still accept
#' R's usual abbreviations.
#'
#' What becomes of a name that matches no argument depends on the interface. The
#' function interface has no estimating equation to forward `...` to, so it
#' requires `...` to be empty and reports an unrecognized name as an error
#' rather than silently ignoring it. The formula interface forwards `...` to
#' `.ee` and matches each name against that function's arguments exactly, naming
#' the argument a refused name was probably meant for. A name that merely
#' abbreviates one is refused too, rather than partially matched to it, since a
#' fit that quietly took a misspelling for `weights` reports different numbers
#' and says nothing about it. An `.ee` that takes `...` of its own accepts any
#' name.
#'
#' @param stacked_equations A formula or a function. When a formula, `data` and
#'   `.ee` must also be provided. When a function, it should take a numeric
#'   vector `theta` and return a p-by-n matrix, whose row names name the
#'   parameters when `init` has none; see [estimate()]. A two-level factor or
#'   character response is converted to a 0/1 indicator against its first level;
#'   an `offset()` term in the formula is passed to `.ee` through its `offset`
#'   argument.
#' @param data A data frame (required when `stacked_equations` is a formula).
#' @param .ee An estimating equation function that accepts `theta`, `X`, and
#'   the response as its third argument, plus optionally additional arguments
#'   (required when `stacked_equations` is a formula). The formula response is
#'   passed positionally, so it reaches whatever the function calls that
#'   argument (`y` for [ee_regression] or [ee_glm], `time` for [ee_aft]). An
#'   equation whose arguments leave any of those nowhere to go cannot be driven
#'   by a formula and is refused before anything is estimated:
#'   [ee_survival_model] takes no design matrix, so it is fitted through the
#'   function interface instead.
#' @param ... For the formula interface, additional arguments passed to `.ee`.
#'   These are evaluated with tidy evaluation in the context of `data`, so
#'   column names can be used directly (e.g., `event = status`). If the model
#'   frame drops rows for missing data, any such argument that spans the full
#'   data is subset to the same rows so it stays aligned with the design matrix
#'   and response. The function interface forwards nothing, so it requires `...`
#'   to be empty.
#' @param init Numeric vector of initial parameter values. When `NULL`
#'   (default) and using the formula interface, a zero vector with names from
#'   the model matrix columns is generated automatically. Names on it label the
#'   parameters and take precedence over the row names of `stacked_equations`.
#'   An explicit `init` with no names of its own takes the same model matrix
#'   names on the formula interface, together with the parameter the estimating
#'   equation estimates beyond the design coefficients (`log_shape` for
#'   [ee_glm] with `"gamma"`, `log_dispersion` with `"negative_binomial"`,
#'   `log_inv_scale` for a non-exponential [ee_aft], `log_sigma` for
#'   [ee_tobit], and `log_phi` for [ee_beta_regression]). An `init` of any
#'   other length is left unnamed, which passes the labeling to the row names
#'   of the estimating functions; see [estimate()] for that channel and for
#'   when the parameters are numbered instead.
#' @param subset Integer vector of parameter indices to solve for, or `NULL`
#'   (default) to solve for all parameters. Indices are 1-based; parameters not
#'   listed are held fixed at their `init` values while the rest are solved.
#'   The variance estimator ignores `subset`. Passed straight through to the
#'   estimator constructor.
#' @param finite_correction Character string for finite-sample correction
#'   (e.g., `"HC1"`), or `NULL` (default) for no correction. When set, the meat
#'   matrix is rescaled and inference switches to the t-distribution with
#'   `df = n_obs - n_params`. Passed straight through to the estimator
#'   constructor.
#' @param solver Character string or function for the solver. Default `NULL`
#'   uses `"rootSolve"` for M-estimation and `"BFGS"` for GMM. See [estimate()]
#'   for the full list of solvers and for how the returned point is judged
#'   against the estimating equations.
#' @param maxiter Integer maximum iterations (default 5000).
#' @param tolerance Numeric convergence tolerance (default 1e-9).
#' @param deriv_method Character string for numerical differentiation method
#'   (default `"capprox"`).
#' @param dx Numeric step size for differentiation (default 1e-9). Must be a
#'   single positive finite number, which is checked whichever `deriv_method` is
#'   in force. The step is absolute and is floored at the floating-point
#'   resolution of each estimate, so a large parameter magnitude cannot silently
#'   reduce it to nothing; see [approx_differentiation()].
#' @param allow_pinv Logical. Use pseudo-inverse if bread is singular? Default
#'   `TRUE`.
#'
#' @returns A fitted `MEstimator` object with populated `theta`, `variance`,
#'   etc. Use [`coef()`][deli-generics], [`vcov()`][deli-generics],
#'   [`confint()`][deli-generics], [`summary()`][deli-display], or
#'   [generics::tidy()] to extract results.
#'
#' @export
#' @examples
#' # Formula interface
#' m <- m_estimate(mpg ~ wt + hp, data = mtcars,
#'                 .ee = ee_regression, model = "linear")
#' coef(m)
#' summary(m)
#'
#' # Function interface
#' y <- c(1, 2, 3, 4, 5)
#' m2 <- m_estimate(
#'   stacked_equations = function(theta) matrix(y - theta[1], nrow = 1),
#'   init = c(mean = 0)
#' )
#' coef(m2)
m_estimate <- function(stacked_equations, ...) {
  UseMethod("m_estimate")
}

#' @rdname m_estimate
#' @export
m_estimate.formula <- function(
  stacked_equations,
  data,
  .ee,
  ...,
  init = NULL,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE
) {
  prep <- prepare_formula_psi(
    stacked_equations,
    data = data,
    .ee = .ee,
    dots = rlang::enquos(...),
    init = init,
    ee_name = formula_ee_name(substitute(.ee)),
    error_call = rlang::current_env()
  )

  obj <- MEstimator(
    stacked_equations = prep$psi,
    init = prep$init,
    subset = subset,
    finite_correction = finite_correction
  )
  obj@model_spec <- prep$model_spec
  estimate(
    obj,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  )
}

#' @rdname m_estimate
#' @export
m_estimate.default <- function(
  stacked_equations,
  ...,
  init,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE
) {
  rlang::check_dots_empty()
  obj <- MEstimator(
    stacked_equations = stacked_equations,
    init = init,
    subset = subset,
    finite_correction = finite_correction
  )
  estimate(
    obj,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  )
}

#' One-step GMM estimation
#'
#' Creates a `GMMEstimator` and estimates it in one call. Supports both a
#' formula interface and a function interface, parallel to [m_estimate()].
#'
#' @details
#' Both interfaces place `...` ahead of `init`, `subset`, `tolerance`, and the
#' rest of the settings, so each of those must be named in full: R does not
#' partially match a supplied name against an argument that follows `...`.
#' [estimate()] and the inference generics take `...` last, so they still accept
#' R's usual abbreviations.
#'
#' What becomes of a name that matches no argument depends on the interface. The
#' function interface has no estimating equation to forward `...` to, so it
#' requires `...` to be empty and reports an unrecognized name as an error
#' rather than silently ignoring it. The formula interface forwards `...` to
#' `.ee` and matches each name against that function's arguments exactly, naming
#' the argument a refused name was probably meant for. A name that merely
#' abbreviates one is refused too, rather than partially matched to it, since a
#' fit that quietly took a misspelling for `weights` reports different numbers
#' and says nothing about it. An `.ee` that takes `...` of its own accepts any
#' name.
#'
#' @inheritParams m_estimate
#' @param overid_maxiter Integer maximum iterations for the two-step iterative
#'   procedure for over-identified problems. Default `200L`. The update converges
#'   linearly rather than quadratically, so a well-identified system commonly
#'   needs tens of passes to reach `overid_tolerance` and a weakly identified one
#'   can need hundreds.
#' @param overid_tolerance Numeric tolerance for convergence of the two-step
#'   iterative procedure. Default `1e-9`.
#'
#' @returns A fitted `GMMEstimator` object.
#'
#' @export
#' @examples
#' # Two instruments for a single treatment effect, confounded by an
#' # unmeasured U. Two moment conditions for one parameter leave the system
#' # over-identified, which is the case GMM is for: `m_estimate()` requires one
#' # estimating equation per parameter.
#' set.seed(42)
#' n <- 200
#' d <- data.frame(Z1 = rbinom(n, 1, 0.5), Z2 = rnorm(n))
#' U <- rnorm(n)
#' d$A <- 0.5 * d$Z1 + 0.3 * d$Z2 + U + rnorm(n)
#' d$Y <- 2 * d$A - U + rnorm(n)
#'
#' # One moment condition per instrument: an instrument should be uncorrelated
#' # with the residual of the outcome on the treatment.
#' ee_iv_moments <- function(theta, X, y, Z) {
#'   t(Z * (y - as.numeric(X %*% theta)))
#' }
#'
#' # The formula interface reads the design and the starting values off the
#' # model and passes anything else, here the instruments, on to `.ee`.
#' g <- gmm_estimate(
#'   Y ~ A - 1,
#'   data = d,
#'   .ee = ee_iv_moments,
#'   Z = cbind(Z1, Z2)
#' )
#'
#' # The instruments move the estimate toward the treatment effect of 2 that
#' # generated the data. The least-squares fit of Y on A ignores the confounding
#' # and stays further from it.
#' coef(g)
#' coef(lm(Y ~ A - 1, data = d))
#'
#' # The function interface takes a `stacked_equations` closure instead, for a
#' # system no formula describes. It has no design to read parameter names from,
#' # so names on `init` are what label the results.
#' psi_iv <- function(theta) {
#'   residual <- d$Y - theta[1] * d$A
#'   rbind(d$Z1 * residual, d$Z2 * residual)
#' }
#'
#' g2 <- gmm_estimate(
#'   stacked_equations = psi_iv,
#'   init = c(effect = 0)
#' )
#' coef(g2)
gmm_estimate <- function(stacked_equations, ...) {
  UseMethod("gmm_estimate")
}

#' @rdname gmm_estimate
#' @export
gmm_estimate.formula <- function(
  stacked_equations,
  data,
  .ee,
  ...,
  init = NULL,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  overid_maxiter = 200L,
  overid_tolerance = 1e-9
) {
  prep <- prepare_formula_psi(
    stacked_equations,
    data = data,
    .ee = .ee,
    dots = rlang::enquos(...),
    init = init,
    ee_name = formula_ee_name(substitute(.ee)),
    error_call = rlang::current_env()
  )

  obj <- GMMEstimator(
    stacked_equations = prep$psi,
    init = prep$init,
    subset = subset,
    finite_correction = finite_correction,
    overid_maxiter = overid_maxiter,
    overid_tolerance = overid_tolerance
  )
  obj@model_spec <- prep$model_spec
  estimate(
    obj,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  )
}

#' @rdname gmm_estimate
#' @export
gmm_estimate.default <- function(
  stacked_equations,
  ...,
  init,
  subset = NULL,
  finite_correction = NULL,
  solver = NULL,
  maxiter = 5000,
  tolerance = 1e-9,
  deriv_method = "capprox",
  dx = 1e-9,
  allow_pinv = TRUE,
  overid_maxiter = 200L,
  overid_tolerance = 1e-9
) {
  rlang::check_dots_empty()
  obj <- GMMEstimator(
    stacked_equations = stacked_equations,
    init = init,
    subset = subset,
    finite_correction = finite_correction,
    overid_maxiter = overid_maxiter,
    overid_tolerance = overid_tolerance
  )
  estimate(
    obj,
    solver = solver,
    maxiter = maxiter,
    tolerance = tolerance,
    deriv_method = deriv_method,
    dx = dx,
    allow_pinv = allow_pinv
  )
}

#' Build the estimating-function closure for the formula interface
#'
#' Shared by [m_estimate.formula()] and [gmm_estimate.formula()]. Constructs the
#' model frame, coerces the response for the regression contract, aligns
#' dots-supplied arguments with the NA-filtered model frame, honors an
#' `offset()` term, checks that `.ee` can be driven by a formula at all and that
#' the names in `...` are arguments it takes, auto-generates `init` when absent,
#' and returns the closure plus the resolved `init`.
#'
#' @param formula The model formula.
#' @param data The data frame.
#' @param .ee The estimating-equation function.
#' @param dots A list of quosures captured from `...`.
#' @param init The user-supplied `init`, or `NULL` to auto-generate it.
#' @param ee_name The name `.ee` was passed under, for the reports about it, or
#'   `NULL` when it arrived as an anonymous function.
#' @param error_call The frame of the method that was called, which every abort
#'   raised here or below reports as the call, since none of the frames they are
#'   raised in belongs to a function the caller wrote.
#'
#' @returns A list with `psi` (the closure), `init` (the resolved vector), and
#'   `model_spec` (what the fit was specified as; see `formula_model_spec()`).
#' @noRd
prepare_formula_psi <- function(
  formula,
  data,
  .ee,
  dots,
  init,
  ee_name = NULL,
  error_call = NULL
) {
  mf <- stats::model.frame(formula, data = data)
  response <- stats::model.response(mf)
  y <- coerce_formula_response(response, error_call = error_call)
  X <- stats::model.matrix(formula, data = mf)

  # Evaluate ... with tidy evaluation in data context.
  ee_args <- lapply(dots, function(q) rlang::eval_tidy(q, data = data))

  # model.frame() may drop rows with missing data. Any dots-supplied argument
  # evaluated against the original data still spans every row, so drop the same
  # rows to keep it aligned with X and y.
  omit <- attr(mf, "na.action")
  if (!is.null(omit)) {
    ee_args <- lapply(ee_args, align_omitted_rows, omit = omit, n = nrow(data))
  }

  # An offset() term reaches the estimating equation through its offset
  # argument; model.matrix() drops the term from the design.
  formula_offset <- stats::model.offset(mf)
  if (!is.null(formula_offset)) {
    if ("offset" %in% names(ee_args)) {
      cli::cli_abort(
        c(
          "An offset was supplied both in the formula and through {.arg ...}.",
          "i" = "Provide the offset in only one place."
        ),
        call = error_call
      )
    }
    ee_args$offset <- formula_offset
  }

  # `do.call()` takes the name of a function as readily as the function itself,
  # so a character `.ee` reaches the estimating equation too. It is resolved to
  # the function it names before anything is checked, so the checks read one
  # function's arguments rather than declining to read a string's, and the rest of
  # the interface recognizes the equation as well. A string is the name to report
  # the equation by.
  if (is.character(.ee) && length(.ee) == 1L) {
    ee_name <- .ee
    .ee <- rlang::try_fetch(
      match.fun(.ee),
      # The lookup is the interface's own step, so its failure is reported the
      # way the checks below are rather than as the `get()` beneath it.
      error = function(cnd) {
        cli::cli_abort(
          c(
            "{.arg .ee} must be a function or the name of one.",
            "x" = "No function named {.val {ee_name}} was found."
          ),
          call = error_call
        )
      }
    )
  }

  # What the interface fills is checked before what the caller named, so an
  # equation a formula cannot drive is reported as that rather than as a
  # misspelling in a call that was never going to work. Both run before the
  # equation is evaluated, so neither fault reaches it.
  check_formula_ee_signature(
    .ee,
    ee_args = ee_args,
    has_formula_offset = !is.null(formula_offset),
    ee_name = ee_name,
    error_call = error_call
  )
  check_formula_ee_dots(
    .ee,
    ee_args = ee_args,
    ee_name = ee_name,
    error_call = error_call
  )

  # The parameter the estimating equation estimates beyond one coefficient per
  # design column, where it is one of the equations that appends any.
  appended <- appended_param_name(.ee, ee_args)

  init_auto <- is.null(init)
  if (init_auto) {
    init <- stats::setNames(rep(0, ncol(X)), colnames(X))
  } else if (is.null(names(init))) {
    # An explicit `init` with no names of its own is labeled from the same
    # column headings, which are otherwise lost: coerce_design() drops them
    # before the estimating equation sees the design, so nothing downstream can
    # recover them. A length neither the design nor the appended parameter
    # accounts for is left alone.
    names(init) <- design_param_names(colnames(X), appended, length(init))
  }

  # The response is passed positionally so it fills the first argument not named
  # theta or X (`y`, `time`, ...).
  psi <- function(theta) {
    do.call(.ee, c(list(theta = theta, X = X), list(y), ee_args))
  }

  if (init_auto) {
    # Record the starting values the user never chose. estimate() evaluates the
    # estimating function at them anyway, and reads this attribute to recognize
    # a failure there as one the automatic length may explain. The appended
    # parameter travels with them so the diagnostic can name what the automatic
    # length leaves out, and the entry point so it reports the call the caller
    # typed from the frame it is raised in, several below this one.
    attr(psi, "deli_auto_init") <- init
    attr(psi, "deli_auto_init_appended") <- appended
    attr(psi, "deli_auto_init_call") <- error_call
  }

  list(
    psi = psi,
    init = init,
    model_spec = formula_model_spec(
      mf = mf,
      X = X,
      y = y,
      .ee = .ee,
      ee_args = ee_args,
      response_levels = formula_response_levels(response)
    )
  )
}

#' Coerce a formula response for the regression contract
#'
#' Character responses become factors, and two-level factors become a 0/1
#' indicator against the first level (the GLM contract). Responses with a
#' different number of levels are rejected, since the formula interface has no
#' way to encode them.
#'
#' @param y The response extracted with [stats::model.response()].
#' @param error_call The frame to report the error against.
#'
#' @returns A numeric or logical response vector.
#' @noRd
coerce_formula_response <- function(y, error_call = NULL) {
  if (is.character(y)) {
    y <- factor(y)
  }
  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      cli::cli_abort(
        c(
          "The formula interface requires a numeric, logical, or two-level
           factor response.",
          "i" = "The response has {nlevels(y)} level{?s}.",
          "i" = "Multinomial models ({.fn ee_mlogit}) take an indicator-matrix
                 response through the function interface."
        ),
        call = error_call
      )
    }
    y <- as.numeric(y != levels(y)[1])
  }
  y
}

#' The levels a formula response was scored against
#'
#' Answers, for a response that `coerce_formula_response()` turns into a 0/1
#' indicator, which category the 0 stands for and which the 1 does. `NULL` for a
#' response that arrives numeric or logical and is scored against nothing.
#'
#' Kept beside `coerce_formula_response()` because the two have to agree on how
#' a character response becomes a factor: the levels reported here are the ones
#' the coercion scores against, and a change to either has to be made in both.
#' The levels are read here rather than returned by the coercion so that the
#' response handed to the estimating equation stays a plain numeric vector,
#' carrying no attribute that could travel into an estimating function's return.
#'
#' @param y The response extracted with [stats::model.response()], before
#'   coercion.
#'
#' @returns A character vector of levels, or `NULL`.
#' @noRd
formula_response_levels <- function(y) {
  if (is.character(y)) {
    y <- factor(y)
  }
  if (is.factor(y)) levels(y) else NULL
}

#' Drop NA-omitted rows from a dots-supplied argument
#'
#' Removes the rows recorded in a model frame's `na.action` attribute from any
#' argument that still spans the original data, leaving arguments of other
#' lengths untouched for downstream validation to catch.
#'
#' @param a A dots-supplied argument (vector or matrix).
#' @param omit The `na.action` attribute of the model frame.
#' @param n The number of rows in the original data.
#'
#' @returns The argument with omitted rows removed when it spanned the data.
#' @noRd
align_omitted_rows <- function(a, omit, n) {
  if (NROW(a) != n) {
    return(a)
  }
  drop <- as.integer(omit)
  if (is.null(dim(a))) {
    a[-drop]
  } else {
    a[-drop, , drop = FALSE]
  }
}

#' The name `.ee` was written as in the call
#'
#' Read from the expression the caller supplied, so a report about the estimating
#' equation can name it and say which of the arguments in the call is the one to
#' change. A function written out in the call itself has no name to report, and
#' neither does one produced by a call other than `::`, so both give `NULL` and
#' are described by the argument they arrived in instead.
#'
#' Takes the expression rather than reading it here, because `substitute()` has
#' to run in the method the caller reached.
#'
#' @param expr The result of `substitute(.ee)` in a `.formula` method.
#'
#' @returns A string, or `NULL`.
#' @noRd
formula_ee_name <- function(expr) {
  if (is.symbol(expr)) {
    name <- as.character(expr)
    # A missing argument substitutes to the empty symbol, which names nothing.
    return(if (nzchar(name)) name else NULL)
  }
  if (is.call(expr) && identical(expr[[1]], quote(`::`))) {
    return(paste(deparse(expr), collapse = ""))
  }
  NULL
}

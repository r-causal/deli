#' Estimating equations for the g-formula (g-computation)
#'
#' Returns a stacked set of estimating equations for the g-formula. When
#' `X0 = NULL`, estimates a single causal mean under the plan encoded by
#' `X1`. When `X0` is provided, estimates the average causal effect
#' (difference between two plans).
#'
#' @param theta Numeric vector. If `X0 = NULL`, length is `1 + p` (causal
#'   mean + regression coefficients). If `X0` is provided, length is
#'   `3 + p` (ACE, mean under X1, mean under X0, regression coefficients).
#' @param y Numeric vector of n observed outcomes.
#' @param X Numeric n-by-p design matrix (observed data).
#' @param X1 Numeric n-by-p design matrix under action plan 1.
#' @param X0 Optional n-by-p design matrix under action plan 0.
#'   Default `NULL`.
#' @param force_continuous Logical. Force linear regression even when `y`
#'   is binary? Default `FALSE`.
#'
#' @returns A matrix of estimating equation contributions. When `X0 = NULL` the
#'   first row is named `causal_mean`; when `X0` is provided the first three
#'   rows are named `ACE`, `E[Y^1]`, and `E[Y^0]`, where 1 and 0 index the two
#'   plans. The outcome model rows are named `X_1` through `X_p` for the
#'   columns of `X`.
#'
#' @examples
#' # A binary treatment, two confounders, and a continuous outcome whose true
#' # average causal effect is 1.5.
#' set.seed(42)
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rbinom(n, 1, 0.4)
#' A <- rbinom(n, 1, inverse_logit(-0.5 + 0.5 * W1 + 0.3 * W2))
#' Y <- 2 + 1.5 * A + W1 - 0.5 * W2 + rnorm(n)
#'
#' X <- cbind(1, A, W1, W2) # Observed design matrix
#' X1 <- cbind(1, 1, W1, W2) # Everyone treated
#' X0 <- cbind(1, 0, W1, W2) # Everyone untreated
#'
#' psi <- function(theta) ee_gformula(theta, y = Y, X = X, X1 = X1, X0 = X0)
#'
#' # theta holds the average causal effect, the mean under treatment, and the
#' # mean under no treatment, followed by the four outcome model coefficients.
#' m <- m_estimate(stacked_equations = psi, init = rep(0, 7))
#' coef(m)[1:3]
#'
#' @export
ee_gformula <- function(theta, y, X, X1, X0 = NULL, force_continuous = FALSE) {
  X <- coerce_design(X)
  y <- as.numeric(y)
  X1 <- coerce_design(X1)
  n <- nrow(X)

  # Reject misaligned plan designs with a clear message rather than the opaque
  # error the downstream rbind would raise, matching Python's up-front checks.
  check_design_dims_match(X, X1, "X", "X1")

  if (!is.null(X0)) {
    X0 <- coerce_design(X0)
    check_design_dims_match(X, X0, "X", "X0")
  }

  # Determine model type based on outcome
  if (all(y %in% c(0, 1)) && !force_continuous) {
    model <- "logistic"
    transform_fn <- inverse_logit
  } else {
    model <- "linear"
    transform_fn <- identity
  }

  if (is.null(X0)) {
    # Single plan: theta = c(mu1, beta)
    mu1 <- theta[1]
    beta <- theta[-1]

    # Outcome regression
    preds_reg <- ee_regression(theta = beta, X = X, y = y, model = model)

    # Mean under X1
    ya1 <- transform_fn(pt_as_vector(X1 %*% beta)) - mu1

    # Stack: (1+p)-by-n
    out <- rbind(matrix(ya1, nrow = 1), preds_reg)
    # The two branches name the same quantity differently. With two plans the
    # superscripts of `E[Y^1]` and `E[Y^0]` index the pair, and the pair is what
    # makes them readable. One plan has no second element to index, so the row
    # takes a bare name for the bare mean it holds.
    rownames(out) <- c("causal_mean", block_param_names("X", ncol(X)))
    out
  } else {
    # Two plans: theta = c(mud, mu1, mu0, beta)
    mud <- theta[1]
    mu1 <- theta[2]
    mu0 <- theta[3]
    beta <- theta[-(1:3)]

    # Outcome regression
    preds_reg <- ee_regression(theta = beta, X = X, y = y, model = model)

    # Mean under X1
    ya1 <- transform_fn(pt_as_vector(X1 %*% beta)) - mu1

    # Mean under X0
    ya0 <- transform_fn(pt_as_vector(X0 %*% beta)) - mu0

    # Average causal effect
    ace <- mu1 - mu0 - mud

    # Stack: (3+p)-by-n
    out <- rbind(
      matrix(ace, nrow = 1, ncol = n),
      matrix(ya1, nrow = 1),
      matrix(ya0, nrow = 1),
      preds_reg
    )
    rownames(out) <- c(
      "ACE",
      "E[Y^1]",
      "E[Y^0]",
      block_param_names("X", ncol(X))
    )
    out
  }
}

#' Estimating equations for inverse probability weighting (IPW)
#'
#' Estimates the average causal effect using IPW with a logistic propensity
#' score model.
#'
#' @param theta Numeric vector of length `3 + b`, where `b` is the number
#'   of propensity score model parameters.
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n binary treatment indicators (0/1).
#' @param W Numeric n-by-b design matrix for the propensity score model.
#' @param truncate Optional length-2 numeric vector `c(lower, upper)` to
#'   clip estimated propensity scores. Bounds must be in ascending order
#'   (`lower <= upper`). Default `NULL`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A `(3+b)`-by-n matrix of estimating equation contributions, with
#'   the first three rows named `ACE`, `E[Y^1]`, and `E[Y^0]` and the propensity
#'   score rows named `W_1` through `W_b` for the columns of `W`.
#'
#' @examples
#' # A binary treatment, two confounders, and a continuous outcome whose true
#' # average causal effect is 1.5.
#' set.seed(42)
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rbinom(n, 1, 0.4)
#' A <- rbinom(n, 1, inverse_logit(-0.5 + 0.5 * W1 + 0.3 * W2))
#' Y <- 2 + 1.5 * A + W1 - 0.5 * W2 + rnorm(n)
#'
#' W_ps <- cbind(1, W1, W2) # Propensity score design matrix
#'
#' psi <- function(theta) ee_ipw(theta, y = Y, A = A, W = W_ps)
#'
#' # theta holds the average causal effect, the mean under treatment, and the
#' # mean under no treatment, followed by the three propensity score
#' # coefficients.
#' m <- m_estimate(stacked_equations = psi, init = rep(0, 6))
#' coef(m)[1:3]
#'
#' @export
ee_ipw <- function(theta, y, A, W, truncate = NULL, weights = NULL) {
  W <- coerce_design(W)
  A <- as.numeric(A)
  y <- as.numeric(y)
  n <- length(y)
  beta <- theta[-(1:3)]

  # Propensity score model (logistic regression of A on W)
  preds_reg <- ee_regression(
    theta = beta,
    X = W,
    y = A,
    model = "logistic"
  )

  # Estimated propensity scores
  pi_hat <- inverse_logit(pt_as_vector(W %*% beta))

  # Truncation
  if (!is.null(truncate)) {
    check_truncate_order(truncate)
    pi_hat <- pmin(pmax(pi_hat, truncate[1]), truncate[2])
  }

  # External weights
  w <- generate_weights(n, weights)

  # Y(a=1): weighted outcome among treated

  ya1 <- (A * y) / pi_hat * w - theta[2]

  # Y(a=0): weighted outcome among untreated
  ya0 <- ((1 - A) * y) / (1 - pi_hat) * w - theta[3]

  # ATE
  ace <- theta[2] - theta[3] - theta[1]

  # Stack: (3+b)-by-n
  out <- rbind(
    matrix(ace, nrow = 1, ncol = n),
    matrix(ya1, nrow = 1),
    matrix(ya0, nrow = 1),
    preds_reg
  )
  rownames(out) <- c(
    "ACE",
    "E[Y^1]",
    "E[Y^0]",
    block_param_names("W", ncol(W))
  )
  out
}

#' Estimating equations for augmented inverse probability weighting (AIPW)
#'
#' Estimates the average causal effect using AIPW, which combines a
#' propensity score model and an outcome model for doubly-robust estimation.
#'
#' @param theta Numeric vector of length `3 + b + c`, where `b` is the
#'   number of propensity score model parameters and `c` is the number of
#'   outcome model parameters.
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n binary treatment indicators (0/1).
#' @param W Numeric n-by-b design matrix for the propensity score model.
#' @param X Numeric n-by-c design matrix for the outcome model.
#' @param X1 Numeric n-by-c design matrix under A=1 for all units.
#' @param X0 Numeric n-by-c design matrix under A=0 for all units.
#' @param truncate Optional length-2 numeric vector `c(lower, upper)` to
#'   clip propensity scores. Bounds must be in ascending order
#'   (`lower <= upper`). Default `NULL`.
#' @param force_continuous Logical. Force linear regression for outcome
#'   model? Default `FALSE`.
#'
#' @returns A `(3+b+c)`-by-n matrix of estimating equation contributions, with
#'   the first three rows named `ACE`, `E[Y^1]`, and `E[Y^0]`, the propensity
#'   score rows named `W_1` through `W_b`, and the outcome model rows named
#'   `X_1` through `X_c`.
#'
#' @examples
#' # A binary treatment, two confounders, and a continuous outcome whose true
#' # average causal effect is 1.5.
#' set.seed(42)
#' n <- 1000
#' W1 <- rnorm(n)
#' W2 <- rbinom(n, 1, 0.4)
#' A <- rbinom(n, 1, inverse_logit(-0.5 + 0.5 * W1 + 0.3 * W2))
#' Y <- 2 + 1.5 * A + W1 - 0.5 * W2 + rnorm(n)
#'
#' W_ps <- cbind(1, W1, W2) # Propensity score design matrix
#' X <- cbind(1, A, W1, W2) # Outcome model, observed design matrix
#' X1 <- cbind(1, 1, W1, W2) # Outcome model, everyone treated
#' X0 <- cbind(1, 0, W1, W2) # Outcome model, everyone untreated
#'
#' psi <- function(theta) {
#'   ee_aipw(theta, y = Y, A = A, W = W_ps, X = X, X1 = X1, X0 = X0)
#' }
#'
#' # theta holds the average causal effect, the mean under treatment, and the
#' # mean under no treatment, followed by the three propensity score
#' # coefficients and the four outcome model coefficients.
#' m <- m_estimate(stacked_equations = psi, init = rep(0, 10))
#' summary(m, subset = 1:3)
#'
#' @export
ee_aipw <- function(
  theta,
  y,
  A,
  W,
  X,
  X1,
  X0,
  truncate = NULL,
  force_continuous = FALSE
) {
  y <- as.numeric(y)
  A <- as.numeric(A)
  W <- coerce_design(W)
  X <- coerce_design(X)
  X1 <- coerce_design(X1)
  X0 <- coerce_design(X0)
  n <- length(y)

  # Misaligned counterfactual designs would otherwise recycle silently into the
  # potential-outcome rows, so reject them up front as Python does.
  check_design_dims_match(X, X1, "X", "X1")
  check_design_dims_match(X, X0, "X", "X0")

  b <- ncol(W)

  # Extract parameters
  mud <- theta[1]
  mu1 <- theta[2]
  mu0 <- theta[3]
  alpha <- theta[4:(3 + b)]
  beta <- theta[(4 + b):length(theta)]

  # Propensity score model
  pi_model <- ee_regression(
    theta = alpha,
    X = W,
    y = A,
    model = "logistic"
  )
  pi_hat <- inverse_logit(pt_as_vector(W %*% alpha))

  # Truncation
  if (!is.null(truncate)) {
    check_truncate_order(truncate)
    pi_hat <- pmin(pmax(pi_hat, truncate[1]), truncate[2])
  }

  # Outcome model
  if (all(y %in% c(0, 1)) && !force_continuous) {
    model <- "logistic"
    transform_fn <- inverse_logit
  } else {
    model <- "linear"
    transform_fn <- identity
  }

  m_model <- ee_regression(theta = beta, y = y, X = X, model = model)
  ya1 <- transform_fn(pt_as_vector(X1 %*% beta))
  ya0 <- transform_fn(pt_as_vector(X0 %*% beta))

  # AIPW estimator
  ace <- mu1 - mu0 - mud
  y1_star <- (y * A / pi_hat - ya1 * (A - pi_hat) / pi_hat) - mu1
  y0_star <- (y * (1 - A) / (1 - pi_hat) + ya0 * (A - pi_hat) / (1 - pi_hat)) -
    mu0

  # Stack: (3+b+c)-by-n
  out <- rbind(
    matrix(ace, nrow = 1, ncol = n),
    matrix(y1_star, nrow = 1),
    matrix(y0_star, nrow = 1),
    pi_model,
    m_model
  )
  rownames(out) <- c(
    "ACE",
    "E[Y^1]",
    "E[Y^0]",
    block_param_names("W", b),
    block_param_names("X", ncol(X))
  )
  out
}

#' Estimating equations for IPW marginal structural model
#'
#' Estimates the parameters of a marginal structural model using inverse
#' probability weighting with a logistic propensity score model.
#'
#' @param theta Numeric vector of length `c + b`, where `c` is the number
#'   of MSM parameters and `b` is the number of propensity score model
#'   parameters.
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n binary treatment indicators (0/1).
#' @param W Numeric n-by-b design matrix for the propensity score model.
#' @param V Numeric n-by-c design matrix for the marginal structural model.
#' @param distribution Character string for the GLM distribution. See
#'   [ee_glm()] for options.
#' @param link Character string for the GLM link function. See [ee_glm()]
#'   for options.
#' @param hyperparameter Optional numeric hyperparameter passed straight
#'   through to the marginal structural model's [ee_glm()] call. Used by the
#'   tweedie outcome distribution, where it fixes the variance power
#'   `v(mu) = mu^hyperparameter`. Default `NULL`. Note that the theta
#'   partition reserves exactly `ncol(V)` slots for the MSM and no slot for
#'   an estimated nuisance parameter, so outcome families that carry one
#'   (`"gamma"`, `"negative_binomial"`, `"nb"`) cannot be used as the MSM
#'   outcome model. Naming one is refused by name, before the outcome model is
#'   formed and whatever `hyperparameter` is set to. Python Delicatessen cannot
#'   fit those families either, where the attempt fails as a shape error.
#' @param truncate Optional length-2 numeric vector `c(lower, upper)` to
#'   clip estimated propensity scores. Bounds must be in ascending order
#'   (`lower <= upper`). Default `NULL`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A `(c+b)`-by-n matrix of estimating equation contributions. The
#'   marginal structural model rows are named `MSM alpha_0` through
#'   `MSM alpha_(c-1)`, matching the zero-based subscripts the literature gives
#'   those parameters. The propensity score rows are named `W_1` through `W_b`
#'   for the columns of `W`.
#'
#' @examples
#' # A confounded binary treatment whose true effect on the outcome is -2.
#' set.seed(42)
#' n <- 500
#' W <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, 0.25 + 0.5 * W)
#' Y <- 5 + 2 * W - 2 * A + rnorm(n)
#'
#' W_ps <- cbind(1, W) # Propensity score design matrix
#' V <- cbind(1, A) # Marginal structural model design matrix
#'
#' psi <- function(theta) {
#'   ee_ipw_msm(
#'     theta,
#'     y = Y,
#'     A = A,
#'     W = W_ps,
#'     V = V,
#'     distribution = "normal",
#'     link = "identity"
#'   )
#' }
#'
#' # theta holds the two marginal structural model coefficients, whose slope is
#' # the causal effect, followed by the two propensity score coefficients.
#' m <- m_estimate(stacked_equations = psi, init = rep(0, 4))
#' coef(m)
#'
#' @export
ee_ipw_msm <- function(
  theta,
  y,
  A,
  W,
  V,
  distribution,
  link,
  hyperparameter = NULL,
  truncate = NULL,
  weights = NULL
) {
  # Coerce inputs
  W <- coerce_design(W)
  V <- coerce_design(V)
  A <- as.numeric(A)
  y <- as.numeric(y)
  n <- length(y)

  # The name is judged before it is compared, because the comparison below reads
  # it through an `if`, which is where a `NULL` or a longer vector failed as
  # base R rather than as an argument of this function.
  check_family_name(distribution, "distribution")
  check_family_name(link, "link")

  # An outcome family that estimates a nuisance parameter has nowhere to keep it
  # here. ee_glm() reads the last element of the coefficient vector it is handed
  # as that parameter, and the theta partition below reserves exactly ncol(V)
  # slots for the MSM, so the design would be multiplied by a coefficient vector
  # one element short. Refusing by name is the point: the partition is this
  # function's, and the failure otherwise surfaces as a non-conformable matrix
  # product the caller never wrote. Lowercased as ee_glm() lowercases it, so a
  # capitalized name is refused here rather than reaching that product.
  if (tolower(distribution) %in% c("gamma", "negative_binomial", "nb")) {
    supported <- c("normal", "poisson", "binomial", "inverse_normal", "tweedie")
    cli::cli_abort(c(
      "The {.val {distribution}} distribution cannot be the outcome model of a
       marginal structural model.",
      "x" = "It estimates a nuisance parameter, and {.arg theta} reserves exactly
             {ncol(V)} slot{?s} for the marginal structural model, one per column
             of {.arg V}, with none for that parameter.",
      "i" = "Supported {.arg distribution} values are {.or {.val {supported}}}."
    ))
  }

  # Extract parameters: first c for MSM, remaining b for PS model
  c_params <- ncol(V)
  alpha <- theta[1:c_params] # MSM parameters
  beta <- theta[(c_params + 1):length(theta)] # PS model parameters

  # Propensity score model (logistic regression of A on W)
  preds_reg <- ee_regression(
    theta = beta,
    X = W,
    y = A,
    model = "logistic",
    weights = NULL
  )

  # Estimated propensity scores
  pi_hat <- inverse_logit(pt_as_vector(W %*% beta))

  # Truncation
  if (!is.null(truncate)) {
    check_truncate_order(truncate)
    pi_hat <- pmin(pmax(pi_hat, truncate[1]), truncate[2])
  }

  # IPW weights: 1/pi for treated, 1/(1-pi) for untreated

  ipw <- ifelse(A == 1, 1 / pi_hat, 1 / (1 - pi_hat))

  # Apply external weights if provided
  if (!is.null(weights)) {
    ipw <- ipw * generate_weights(n, weights)
  }

  # Marginal structural model via weighted GLM
  ee_msm <- ee_glm(
    theta = alpha,
    X = V,
    y = y,
    distribution = distribution,
    link = link,
    hyperparameter = hyperparameter, # Passed through for tweedie variance power
    weights = ipw,
    offset = NULL
  )

  # Stack: (c+b)-by-n
  out <- rbind(ee_msm, preds_reg)
  # The marginal structural model coefficients are subscripted from zero, as
  # alpha_0 and alpha_1 are in the literature this implements, so the sequence
  # starts one below the row position. The `hyperparameter` documentation
  # explains why the MSM block is always exactly ncol(V) rows: the outcome
  # families that carry an estimated nuisance parameter have no slot in the
  # theta partition and cannot reach this stack.
  rownames(out) <- c(
    sprintf("MSM alpha_%d", seq_len(c_params) - 1L),
    block_param_names("W", ncol(W))
  )
  out
}

#' Estimating equations for g-estimation of structural nested mean models
#'
#' Estimates the parameters of a structural nested mean model via
#' g-estimation. Supports both inefficient (X = NULL) and efficient
#' (X provided) g-estimators, and both linear and log-linear (Poisson)
#' structural mean models.
#'
#' @param theta Numeric vector. For the inefficient g-estimator,
#'   length is `b + c` (SMM parameters + PS model parameters). For the
#'   efficient g-estimator, length is `b + c + d` (SMM + PS + outcome
#'   model parameters).
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n binary treatment indicators (0/1).
#' @param W Numeric n-by-c design matrix for the propensity score model.
#' @param V Numeric n-by-b design matrix for the structural mean model.
#'   Should NOT include A itself.
#' @param X Optional n-by-d design matrix for the outcome model
#'   (efficient g-estimator). Default `NULL` (inefficient g-estimator).
#' @param model Character string: `"linear"` or `"poisson"`. Default
#'   `"linear"`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A matrix of estimating equation contributions. The structural mean
#'   model rows are named `SNM phi_0` through `SNM phi_(b-1)`, matching the
#'   zero-based subscripts the literature gives those parameters. The propensity
#'   score rows are named `W_1` through `W_c`, and, for the efficient
#'   g-estimator, the outcome model rows are named `X_1` through `X_d`.
#'
#' @examples
#' # A confounded binary treatment whose true effect on the outcome is -2.
#' set.seed(42)
#' n <- 500
#' W <- rbinom(n, 1, 0.5)
#' A <- rbinom(n, 1, 0.25 + 0.5 * W)
#' Y <- 5 + 2 * W - 2 * A + rnorm(n)
#'
#' W_ps <- cbind(1, W) # Propensity score design matrix
#'
#' # An intercept-only structural mean model gives a single causal contrast.
#' # Build it with rep() so the column has one entry per observation.
#' V <- cbind(rep(1, n))
#'
#' psi <- function(theta) {
#'   ee_gestimation_snmm(theta, y = Y, A = A, W = W_ps, V = V, model = "linear")
#' }
#'
#' # theta holds the structural mean model coefficient, which is the causal
#' # effect, followed by the two propensity score coefficients.
#' m <- m_estimate(stacked_equations = psi, init = rep(0, 3))
#' coef(m)
#'
#' @export
ee_gestimation_snmm <- function(
  theta,
  y,
  A,
  W,
  V,
  X = NULL,
  model = "linear",
  weights = NULL
) {
  # Coerce inputs
  y <- as.numeric(y)
  A <- as.numeric(A)
  W <- coerce_design(W)
  V <- coerce_design(V)
  n <- length(y)

  # Parameter indexing
  pdiv <- ncol(V) # Number of SMM parameters
  qdiv <- pdiv + ncol(W) # End of PS model parameters

  # Process weights
  w <- generate_weights(n, weights)

  # Extract parameter subsets
  phi <- theta[1:pdiv] # SMM parameters
  alpha <- theta[(pdiv + 1):qdiv] # PS model parameters
  if (!is.null(X)) {
    X <- coerce_design(X)
    beta <- theta[(qdiv + 1):length(theta)] # Outcome model parameters
  }

  # Compute H(phi) based on model type
  model <- tolower(model)
  # V * A: element-wise multiplication of each column of V by A
  VA <- V * A
  if (model == "linear") {
    # H(phi) = Y - (V*A) %*% phi
    h_phi <- y - pt_as_vector(VA %*% phi)
    y_transform <- identity_transform
  } else if (model == "poisson") {
    # H(phi) = Y * exp(-(V*A) %*% phi)
    h_phi <- y * exp(-pt_as_vector(VA %*% phi))
    y_transform <- exp
  } else {
    cli::cli_abort(
      "model={.val {model}} is not supported. Use {.val linear} or {.val poisson}."
    )
  }

  # Propensity score model: E[A | W]
  ee_log <- ee_regression(
    theta = alpha,
    X = W,
    y = A,
    model = "logistic",
    weights = weights
  )
  pi_hat <- inverse_logit(pt_as_vector(W %*% alpha))
  a_resid <- A - pi_hat # Residuals for A

  # Estimating equations for the g-estimator

  if (!is.null(X)) {
    # Efficient g-estimator: include outcome model E[H(phi) | W]
    ee_out <- ee_regression(
      theta = beta,
      X = X,
      y = h_phi,
      model = model,
      weights = weights
    )
    yhat <- y_transform(pt_as_vector(X %*% beta)) # Predicted H(phi)
  } else {
    # Inefficient g-estimator: no outcome model
    yhat <- 0
  }

  # Structural mean model estimating equation
  y0_resid <- h_phi - yhat
  # ee_smm: b-by-n matrix
  ee_smm <- t(V * (w * a_resid * y0_resid))

  # Stack estimating equations. The structural mean model parameters are
  # subscripted from zero, as phi_0 and phi_1 are in the literature this
  # implements, so the sequence starts one below the row position.
  smm_names <- sprintf("SNM phi_%d", seq_len(pdiv) - 1L)
  ps_names <- block_param_names("W", ncol(W))
  if (!is.null(X)) {
    out <- rbind(ee_smm, ee_log, ee_out)
    rownames(out) <- c(smm_names, ps_names, block_param_names("X", ncol(X)))
  } else {
    out <- rbind(ee_smm, ee_log)
    rownames(out) <- c(smm_names, ps_names)
  }
  out
}

#' Estimating equations for instrumental variable (IV) estimation
#'
#' Estimates the causal effect using the usual IV / Wald estimator.
#' The parameter of interest is the additive effect of treatment A on
#' outcome Y leveraging instrument Z.
#'
#' @param theta Numeric vector of length 2: the causal effect and the
#'   mean of the instrument.
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n observed treatment values.
#' @param Z Numeric vector of n binary instrument values (0/1).
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A 2-by-n matrix of estimating equation contributions, with rows
#'   named `causal_effect` and `mean_Z`.
#'
#' @examples
#' # An unmeasured confounder U biases the association between A and Y, but the
#' # instrument Z affects Y only through A. The true causal effect is 3.
#' set.seed(123)
#' n <- 500
#' Z <- rbinom(n, 1, 0.5)
#' U <- rnorm(n)
#' A <- rbinom(n, 1, inverse_logit(-1 + 3 * Z + U))
#' Y <- 3 * A - U + rnorm(n, sd = 0.5)
#'
#' psi <- function(theta) ee_iv_causal(theta, y = Y, A = A, Z = Z)
#'
#' # theta holds the causal effect followed by the mean of the instrument.
#' m <- m_estimate(stacked_equations = psi, init = c(0, 0.5))
#' coef(m)
#'
#' @export
ee_iv_causal <- function(theta, y, A, Z, weights = NULL) {
  # Coerce inputs
  y <- as.numeric(y)
  A <- as.numeric(A)
  Z <- as.numeric(Z)
  n <- length(y)

  # Process weights
  w <- generate_weights(n, weights)

  # Usual IV estimating equations
  # theta[1] = causal effect (beta)
  # theta[2] = mean of Z (mu)
  ee_prz <- w * (Z - theta[2])
  ee_iva <- w * (y - theta[1] * A) * (Z - theta[2])

  # Stack: 2-by-n
  out <- rbind(
    matrix(ee_iva, nrow = 1),
    matrix(ee_prz, nrow = 1)
  )
  rownames(out) <- c("causal_effect", "mean_Z")
  out
}

#' Estimating equations for Two-Stage Least Squares (2SLS)
#'
#' Estimates the causal effect using two-stage least squares for
#' instrumental variable analysis. The first stage regresses the
#' treatment on the instruments (and exogenous variables), and the
#' second stage regresses the outcome on predicted treatment (and
#' exogenous variables).
#'
#' @param theta Numeric vector of length `1 + b + 2c`, where `b` is
#'   the number of instruments in `Z` and `c` is the number of
#'   exogenous variables in `W`. The first `1 + c` parameters are for
#'   the second-stage model; the remainder are for the first-stage model.
#' @param y Numeric vector of n observed outcomes.
#' @param A Numeric vector of n observed treatment values.
#' @param Z Numeric n-by-b matrix of instrumental variable(s).
#' @param W Optional n-by-c matrix of exogenous variables included in
#'   both stages. Default `NULL`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#'
#' @returns A `(1+b+2c)`-by-n matrix of estimating equation contributions. The
#'   second-stage rows are named `stage2_A` for the fitted treatment and
#'   `stage2_W_1` through `stage2_W_c`; the first-stage rows are named
#'   `stage1_Z_1` through `stage1_Z_b` and `stage1_W_1` through `stage1_W_c`.
#'   The columns of `W` appear in both stages, so each name carries the stage it
#'   belongs to.
#'
#' @examples
#' # A continuous treatment confounded by an unmeasured U, a strong instrument
#' # Z, and one measured exogenous covariate. The true causal effect is 2.
#' set.seed(42)
#' n <- 500
#' W1 <- rnorm(n)
#' Z <- cbind(rbinom(n, 1, 0.5))
#' U <- rnorm(n)
#' A <- 1.5 * Z[, 1] + 0.3 * W1 + U + rnorm(n, sd = 0.5)
#' Y <- 2 * A - U + 0.5 * W1 + rnorm(n)
#' W <- cbind(1, W1)
#'
#' # At a first stage of all zeros the fitted treatment is identically zero, so
#' # the second stage degenerates: its leading design column vanishes and the
#' # coefficient on it, the causal effect, has nothing left to be estimated
#' # from. Seed the starting values with the ordinary least squares fit of each
#' # stage instead: the first stage regresses A on the instrument and the
#' # exogenous covariates, the second regresses Y on the fitted A and the same
#' # covariates.
#' alpha_init <- as.numeric(coef(lm(A ~ cbind(Z, W) - 1)))
#' a_hat <- as.numeric(cbind(Z, W) %*% alpha_init)
#' beta_init <- as.numeric(coef(lm(Y ~ cbind(a_hat, W) - 1)))
#'
#' psi <- function(theta) ee_2sls(theta, y = Y, A = A, Z = Z, W = W)
#'
#' # theta holds the three second-stage coefficients followed by the three
#' # first-stage coefficients.
#' m <- m_estimate(
#'   stacked_equations = psi,
#'   init = c(beta_init, alpha_init)
#' )
#'
#' # theta_1 is the coefficient on the fitted treatment, the causal effect.
#' coef(m)
#'
#' @export
ee_2sls <- function(theta, y, A, Z, W = NULL, weights = NULL) {
  # Coerce inputs
  y <- as.numeric(y)
  a <- as.numeric(A)
  Z <- coerce_design(Z)

  # Processing parameter vector
  if (is.null(W)) {
    id2s <- 1 # Split point: no exogenous covariates
  } else {
    W <- coerce_design(W)
    id2s <- 1 + ncol(W) # Split point for first/second stage
  }
  beta <- theta[1:id2s] # Second-stage parameters
  alpha <- theta[(id2s + 1):length(theta)] # First-stage parameters

  # Processing design matrices
  if (!is.null(W)) {
    dmatrix1 <- cbind(Z, W) # Stack instruments and exogenous vars
  } else {
    dmatrix1 <- Z
  }
  a_hat <- pt_as_vector(dmatrix1 %*% alpha) # Predicted values of A

  if (!is.null(W)) {
    dmatrix2 <- cbind(a_hat, W) # Stack predicted A and exogenous vars
  } else {
    dmatrix2 <- matrix(a_hat, ncol = 1)
  }

  # First-stage least squares: A ~ Z (+ W)
  ee_stageone <- ee_regression(
    theta = alpha,
    X = dmatrix1,
    y = a,
    model = "linear",
    weights = weights
  )

  # Second-stage least squares: Y ~ A_hat (+ W)
  ee_stagetwo <- ee_regression(
    theta = beta,
    X = dmatrix2,
    y = y,
    model = "linear",
    weights = weights
  )

  # Output: stack second stage on top of first stage

  out <- rbind(ee_stagetwo, ee_stageone)
  # Every column of W is fitted twice, once per stage, so the stage is part of
  # each name. Without it the two blocks would repeat every W label, and a
  # repeated label costs the whole return its names.
  if (is.null(W)) {
    w2_names <- character(0)
    w1_names <- character(0)
  } else {
    w2_names <- block_param_names("stage2_W", ncol(W))
    w1_names <- block_param_names("stage1_W", ncol(W))
  }
  rownames(out) <- c(
    "stage2_A",
    w2_names,
    block_param_names("stage1_Z", ncol(Z)),
    w1_names
  )
  out
}

#' Estimating equations for weighted sensitivity analysis of the mean
#'
#' Estimates the mean of an outcome with missing data using a weighted
#' sensitivity analysis approach. Handles MCAR, MAR, and MNAR
#' mechanisms by specifying a user-defined sensitivity function `q_eval`
#' and a monotone increasing distribution function `H_function`.
#'
#' @param theta Numeric vector of length `1 + b`, where `b` is the
#'   number of columns in `X`. The first element is the corrected
#'   mean; the remainder are regression coefficients.
#' @param y Numeric vector of n outcome values. Missing values should
#'   be indicated via the `delta` parameter.
#' @param delta Numeric vector of n indicators: 1 if `y` is observed,
#'   0 if missing. Must not contain `NA`.
#' @param X Numeric n-by-b design matrix for the missingness model.
#'   Should include an intercept column. Must not contain `NA`.
#' @param q_eval Numeric vector of n evaluated sensitivity function
#'   values, i.e. `q(Y, alpha)`.
#' @param H_function A function mapping real values to `[0, 1]` that
#'   is monotone increasing (e.g., [inverse_logit()]).
#'
#' @returns A `(1+b)`-by-n matrix of estimating equation contributions, with the
#'   first row named `corrected_mean` and the missingness model rows named `X_1`
#'   through `X_b` for the columns of `X`.
#'
#' @examples
#' # An outcome observed for only part of the sample, with missingness driven by
#' # the measured covariate W.
#' set.seed(42)
#' n <- 500
#' W <- rbinom(n, 1, 0.5)
#' Y_full <- 200 - 35 * W + rnorm(n, sd = 5)
#' delta <- rbinom(n, 1, inverse_logit(2 + W))
#'
#' # Missing outcomes never enter the estimating equation, so any placeholder
#' # value works; a zero keeps the arithmetic finite.
#' Y <- ifelse(delta == 1, Y_full, 0)
#' X <- cbind(1, W) # Missingness model design matrix
#'
#' # A sensitivity function of zero everywhere assumes the outcome is missing at
#' # random given W. Nonzero values encode departures from that assumption.
#' psi <- function(theta) {
#'   ee_mean_sensitivity_analysis(
#'     theta,
#'     y = Y,
#'     delta = delta,
#'     X = X,
#'     q_eval = rep(0, n),
#'     H_function = inverse_logit
#'   )
#' }
#'
#' # theta holds the corrected mean, started near the observed outcomes,
#' # followed by the two missingness model coefficients.
#' m <- m_estimate(stacked_equations = psi, init = c(180, 0, 0))
#' coef(m)
#'
#' @export
ee_mean_sensitivity_analysis <- function(
  theta,
  y,
  delta,
  X,
  q_eval,
  H_function
) {
  # Coerce inputs
  delta <- as.numeric(delta)
  y <- as.numeric(y)
  X <- coerce_design(X)
  qy <- as.numeric(q_eval)
  beta <- theta[-1] # Nuisance parameters

  # Read the column count before the scalar-design branch below, which strips
  # the dim off a one-column design and would leave ncol() answering NULL.
  b <- ncol(X)
  n <- length(y)

  # q_eval and delta are per-observation data, never tangent containers, so a
  # short vector would recycle into the mean and nuisance rows. Reject it as
  # Python's broadcast check does.
  check_data_length(qy, n, "q_eval")
  check_data_length(delta, n, "delta")

  # Predicted values from design matrix and nuisance coefficients
  pred_values <- pt_as_vector(X %*% beta) # Linear predictor

  # Solving for the sensitivity analysis mean
  numerator <- delta * y # Numerator for mean EE
  denominator <- H_function(pred_values + qy) # Denominator for mean EE
  # Set missing Y contributions to zero
  ym_ind <- ifelse(delta == 1, numerator / denominator, 0)
  ef_mean <- ym_ind - theta[1] # Sensitivity analysis EE

  # Solving for intercept and coefficients of model
  h_factor <- delta / denominator - 1 # Length-n weighting residual
  # A scalar intercept design reaches here as a length-1 array. Multiplying the
  # length-n factor by that array triggers R's array-recycling deprecation
  # warning, so strip its dim to make this scalar-vector arithmetic. Tangent
  # containers keep their structure; only a plain length-1 array is coerced.
  if (!is_tangent_container(X) && length(X) == 1L) {
    X <- as.numeric(X)
  }
  ef_H <- h_factor * X

  # Returning stacked estimating equations
  out <- rbind(
    matrix(ef_mean, nrow = 1), # theta[1]: sensitivity mean
    t(ef_H) # theta[2:]: nuisance parameters
  )
  rownames(out) <- c("corrected_mean", block_param_names("X", b))
  out
}

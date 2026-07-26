#' Estimating equation for generalized linear models
#'
#' Returns a p-by-n matrix of estimating equation contributions for a GLM
#' with specified distribution and link function:
#' \deqn{\psi_i(\theta) = \{Y_i - g^{-1}(X_i^T \theta)\} \frac{D(\theta)}{v(\theta)} X_i}
#'
#' @param theta Numeric vector of length p.
#' @param X Numeric n-by-p design matrix.
#' @param y Numeric vector of n observed outcome values.
#' @param distribution Character string: `"normal"` (or `"gaussian"`),
#'   `"binomial"` (or `"bernoulli"`), `"poisson"`, `"gamma"`,
#'   `"negative_binomial"` (or `"nb"`), `"inverse_normal"` (or
#'   `"inverse_gaussian"`), `"tweedie"`.
#' @param link Character string: `"identity"`, `"log"`, `"logit"` (or
#'   `"logistic"`), `"probit"`, `"cauchit"`, `"loglog"`, `"cloglog"`,
#'   `"inverse"`, `"sqrt"` (or `"square_root"`).
#' @param hyperparameter Numeric scalar power \eqn{p} for the variance function
#'   of the `"tweedie"` distribution. Must be non-negative. Ignored by every
#'   other distribution. Default `NULL`.
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @section Gamma:
#' The gamma distribution estimates one additional shape parameter beyond the
#' regression coefficients, so `theta` has length `ncol(X) + 1`. The final
#' element is the log of the shape \eqn{\alpha} (the reciprocal of the GLM
#' dispersion \eqn{\phi = 1 / \alpha}), and the variance function is
#' \eqn{\mu^2}. An extra nuisance estimating equation for \eqn{\log(\alpha)} is
#' appended:
#' \deqn{
#'   \left(1 - \frac{Y_i}{\mu_i}\right)
#'   + \log\left(\frac{\alpha Y_i}{\mu_i}\right)
#'   - \psi(\alpha)
#' }
#' where \eqn{\psi} is the digamma function. The returned matrix therefore
#' has `ncol(X) + 1` rows.
#'
#' @section Negative binomial:
#' The negative binomial distribution estimates one additional dispersion
#' parameter beyond the regression coefficients, so `theta` has length
#' `ncol(X) + 1`. The final element is the log of the dispersion
#' \eqn{\alpha}, and the variance function is \eqn{\mu + \alpha \mu^2}. An
#' extra nuisance estimating equation for \eqn{\log(\alpha)} is appended so
#' that the uncertainty in the dispersion is carried honestly through the
#' sandwich variance:
#' \deqn{
#'   -\alpha^{-2}\psi\left(Y_i + \alpha^{-1}\right)
#'   + \alpha^{-2}\psi\left(\alpha^{-1}\right)
#'   + \frac{Y_i}{\alpha^2 \mu_i + \alpha}
#'   - \frac{\frac{\alpha \mu_i}{\alpha \mu_i + 1}
#'       + \log\left(\frac{1}{\alpha \mu_i + 1}\right)}{\alpha^2}
#' }
#' where \eqn{\psi} is the digamma function. The returned matrix therefore
#' has `ncol(X) + 1` rows.
#'
#' @section Tweedie:
#' The tweedie distribution uses the variance function \eqn{v(\mu) = \mu^p},
#' where the power \eqn{p} is the fixed `hyperparameter` rather than an
#' estimated parameter. No nuisance estimating equation is appended, so `theta`
#' has length `ncol(X)` and the returned matrix has `ncol(X)` rows. The power
#' recovers familiar special cases: \eqn{p = 1} gives the Poisson variance
#' \eqn{\mu} and \eqn{p = 2} gives the gamma variance \eqn{\mu^2}, so those
#' choices share the corresponding beta score equations. The `hyperparameter`
#' must be non-negative.
#'
#' @returns A p-by-n matrix. For the gamma and negative binomial
#'   distributions the matrix has `ncol(X) + 1` rows, the final row being the
#'   nuisance equation for that distribution's extra parameter (the gamma shape
#'   or the negative binomial dispersion).
#'
#' @export
#' @examples
#' # Negative binomial GLM with a log link estimates the regression
#' # coefficients plus one log-dispersion parameter, so init has ncol(X) + 1
#' # entries.
#' n <- 200
#' X <- cbind(1, rnorm(n), rnorm(n))
#' mu <- exp(0.5 + 0.5 * X[, 2] - 0.3 * X[, 3])
#' y <- rnbinom(n, size = 2, mu = mu)
#'
#' psi <- function(theta) {
#'   ee_glm(theta, X = X, y = y, distribution = "negative_binomial",
#'          link = "log")
#' }
#'
#' m <- MEstimator(stacked_equations = psi, init = c(0, 0, 0, 0)) |>
#'   estimate()
#' m@theta
#'
#' # Tweedie GLM with a log link. The power p is a fixed hyperparameter, not an
#' # estimated parameter, so init has ncol(X) entries. Here p = 1.5 sits in the
#' # compound Poisson-gamma regime, appropriate for non-negative data with a
#' # point mass at zero.
#' mu <- exp(0.5 + 0.5 * X[, 2] - 0.3 * X[, 3])
#' y_tw <- rpois(n, lambda = mu)
#'
#' psi_tw <- function(theta) {
#'   ee_glm(theta, X = X, y = y_tw, distribution = "tweedie", link = "log",
#'          hyperparameter = 1.5)
#' }
#'
#' m_tw <- MEstimator(stacked_equations = psi_tw, init = c(0, 0, 0)) |>
#'   estimate()
#' m_tw@theta
ee_glm <- function(
  theta,
  X,
  y,
  distribution,
  link,
  hyperparameter = NULL,
  weights = NULL,
  offset = NULL
) {
  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  w <- generate_weights(n, weights)

  dist <- tolower(distribution)

  # Distributions with an extra parameter carry the log of that parameter as
  # the final element of theta: the gamma shape or the negative binomial
  # dispersion.
  if (dist %in% c("gamma", "negative_binomial", "nb")) {
    beta <- theta[-length(theta)] # Regression coefficients
    alpha <- exp(theta[length(theta)]) # Extra parameter on the natural scale
  } else {
    beta <- theta
    alpha <- NULL
  }

  # Linear predictor
  eta <- pt_as_vector(X %*% beta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }

  # Inverse link and its derivative
  inv <- inverse_link(eta, link)
  pred_y <- inv$mu
  deriv <- inv$dmu

  # Distribution variance
  v <- distribution_variance(
    pred_y,
    distribution,
    alpha = alpha,
    hyperparameter = hyperparameter
  )

  # Score: (y - mu) * d(mu)/d(eta) / V(mu) * X
  residuals <- y - pred_y
  ee_beta <- t(X * (w * residuals * deriv / v))

  if (dist == "gamma") {
    # Nuisance estimating equation for the gamma shape parameter. Unlike the
    # negative binomial nuisance row, this one carries the weights to match
    # Python Delicatessen.
    ee_alpha <- w *
      ((1 - y / pred_y) +
        log(alpha * y / pred_y) -
        deli_polygamma(0, alpha))

    # Stack the shape nuisance row beneath the beta score rows.
    out <- rbind(ee_beta, ee_alpha)
    rownames(out) <- NULL
    return(out)
  }

  if (dist %in% c("negative_binomial", "nb")) {
    # Nuisance estimating equation for the log-dispersion parameter, broken
    # into simpler pieces. This row is built without weights to match Python
    # Delicatessen, which weights the beta scores but not this nuisance row.
    p1 <- -alpha^-2 * deli_polygamma(0, y + 1 / alpha)
    p2 <- alpha^-2 * deli_polygamma(0, 1 / alpha)
    p3 <- y / (alpha^2 * pred_y + alpha)
    p4 <- -(alpha *
      pred_y /
      (alpha * pred_y + 1) +
      log(1 / (alpha * pred_y + 1))) /
      alpha^2
    ee_alpha <- p1 + p2 + p3 + p4

    # Stack the dispersion nuisance row beneath the beta score rows.
    out <- rbind(ee_beta, ee_alpha)
    rownames(out) <- NULL
    return(out)
  }

  ee_beta
}

#' Estimating equation for robust regression
#'
#' Returns a p-by-n matrix for robust regression using the specified loss
#' function (linear regression only):
#' \deqn{\psi_i(\theta) = f_k(Y_i - X_i^T \theta) X_i}
#'
#' @param theta Numeric vector of length p.
#' @param X Numeric n-by-p design matrix.
#' @param y Numeric vector of n observed outcome values.
#' @param model Character string. Currently only `"linear"` is supported.
#' @param k Numeric tuning parameter for the loss function.
#' @param loss Character string specifying the loss function. Default
#'   `"huber"`. See [robust_loss_functions()].
#' @param weights Optional numeric vector of n weights. Default `NULL`.
#' @param offset Optional numeric vector of n offsets. Default `NULL`.
#'
#' @returns A p-by-n matrix.
#'
#' @export
ee_robust_regression <- function(
  theta,
  X,
  y,
  model,
  k,
  loss = "huber",
  weights = NULL,
  offset = NULL
) {
  model <- tolower(model)
  if (model != "linear") {
    cli::cli_abort(
      "Only {.val linear} regression is supported for robust regression."
    )
  }

  X <- coerce_design(X)
  y <- coerce_outcome(y)
  n <- nrow(X)
  w <- generate_weights(n, weights)

  # Linear predictor
  eta <- pt_as_vector(X %*% theta)
  if (!is.null(offset)) {
    eta <- eta + as.numeric(offset)
  }

  # Robust residuals
  residuals <- robust_loss_functions(y - eta, loss = loss, k = k)

  # Score
  t(X * (w * residuals))
}

#' Internal inverse link function dispatcher
#' @noRd
inverse_link <- function(eta, link) {
  link <- tolower(link)
  if (link == "identity") {
    list(mu = eta, dmu = rep(1, length(eta)))
  } else if (link == "log") {
    mu <- exp(eta)
    list(mu = mu, dmu = mu)
  } else if (link %in% c("logit", "logistic")) {
    mu <- 1 / (1 + exp(-eta))
    list(mu = mu, dmu = mu * (1 - mu))
  } else if (link == "probit") {
    list(mu = standard_normal_cdf(eta), dmu = standard_normal_pdf(eta))
  } else if (link %in% c("cauchit", "cauchy")) {
    mu <- (1 / pi) * atan(eta) + 0.5
    dmu <- 1 / (pi * (1 + eta^2))
    list(mu = mu, dmu = dmu)
  } else if (link == "loglog") {
    mu <- exp(-exp(-eta))
    dmu <- -exp(-eta - exp(-eta))
    list(mu = mu, dmu = dmu)
  } else if (link == "cloglog") {
    mu <- 1 - exp(-exp(eta))
    dmu <- exp(eta - exp(eta))
    list(mu = mu, dmu = dmu)
  } else if (link == "inverse") {
    mu <- 1 / eta
    dmu <- -1 / (eta^2)
    list(mu = mu, dmu = dmu)
  } else if (link %in% c("sqrt", "square_root")) {
    mu <- eta^2
    dmu <- 2 * eta
    list(mu = mu, dmu = dmu)
  } else {
    cli::cli_abort("Link function {.val {link}} is not supported.")
  }
}

#' Internal distribution variance function
#' @noRd
distribution_variance <- function(
  mu,
  distribution,
  alpha = NULL,
  hyperparameter = NULL
) {
  dist <- tolower(distribution)
  if (dist %in% c("normal", "gaussian")) {
    rep(1, length(mu))
  } else if (dist == "poisson") {
    mu
  } else if (dist %in% c("binomial", "bin", "bernoulli")) {
    mu * (1 - mu)
  } else if (dist == "gamma") {
    mu^2
  } else if (dist %in% c("negative_binomial", "nb")) {
    mu + alpha * mu^2
  } else if (dist %in% c("inverse_normal", "inverse_gaussian")) {
    mu^3
  } else if (dist == "tweedie") {
    # The tweedie variance is mu raised to the power hyperparameter, which must
    # be supplied and non-negative.
    if (is.null(hyperparameter)) {
      cli::cli_abort(
        "The {.val tweedie} distribution requires a {.arg hyperparameter}."
      )
    }
    if (hyperparameter < 0) {
      cli::cli_abort(
        "The {.arg hyperparameter} for the {.val tweedie} distribution must be non-negative."
      )
    }
    mu^hyperparameter
  } else {
    cli::cli_abort("Distribution {.val {distribution}} is not supported.")
  }
}

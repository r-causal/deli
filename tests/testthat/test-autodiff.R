# Tests for auto_differentiation (forward-mode autodiff)

# ---- scalar function, single input ----

test_that("autodiff computes derivative of x^2", {
  f <- function(x) x[1]^2
  result <- auto_differentiation(c(3), f)
  # d/dx x^2 = 2x = 6
  expect_equal(result[1, 1], 6)
})

test_that("autodiff computes derivative of x^3 - 2x", {
  f <- function(x) x[1]^3 - 2 * x[1]
  result <- auto_differentiation(c(2), f)
  # d/dx(x^3 - 2x) = 3x^2 - 2 = 10
  expect_equal(result[1, 1], 10)
})

test_that("autodiff computes derivative of constant", {
  f <- function(x) 5
  result <- auto_differentiation(c(1), f)
  expect_equal(result[1, 1], 0)
})

# ---- multiple inputs, single output ----

test_that("autodiff partial derivatives of x1*x2", {
  f <- function(x) x[1] * x[2]
  result <- auto_differentiation(c(3, 4), f)
  # df/dx1 = x2 = 4, df/dx2 = x1 = 3
  expect_equal(result[1, 1], 4)
  expect_equal(result[1, 2], 3)
})

test_that("autodiff gradient of x1^2 + x2^2", {
  f <- function(x) x[1]^2 + x[2]^2
  result <- auto_differentiation(c(3, 4), f)
  # df/dx1 = 2*x1 = 6, df/dx2 = 2*x2 = 8
  expect_equal(result[1, 1], 6)
  expect_equal(result[1, 2], 8)
})

# ---- multiple inputs, multiple outputs (Jacobian) ----

test_that("autodiff Jacobian of vector function", {
  f <- function(x) c(x[1]^2 - x[2], x[1] * x[2])
  result <- auto_differentiation(c(3, 4), f)
  # Row 1: d(x1^2 - x2)/dx1 = 2*x1 = 6, d(x1^2 - x2)/dx2 = -1
  # Row 2: d(x1*x2)/dx1 = x2 = 4, d(x1*x2)/dx2 = x1 = 3
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 2)
  expect_equal(result[1, 1], 6)
  expect_equal(result[1, 2], -1)
  expect_equal(result[2, 1], 4)
  expect_equal(result[2, 2], 3)
})

# ---- exp, log, sqrt ----

test_that("autodiff derivative of exp(x)", {
  f <- function(x) exp(x[1])
  result <- auto_differentiation(c(1), f)
  # d/dx exp(x) = exp(x) = exp(1)
  expect_equal(result[1, 1], exp(1))
})

test_that("autodiff derivative of log(x)", {
  f <- function(x) log(x[1])
  result <- auto_differentiation(c(2), f)
  # d/dx log(x) = 1/x = 0.5
  expect_equal(result[1, 1], 0.5)
})

test_that("autodiff derivative of sqrt(x)", {
  f <- function(x) sqrt(x[1])
  result <- auto_differentiation(c(4), f)
  # d/dx sqrt(x) = 1/(2*sqrt(x)) = 0.25
  expect_equal(result[1, 1], 0.25)
})

# ---- trigonometric functions ----

test_that("autodiff derivative of sin(x)", {
  f <- function(x) sin(x[1])
  result <- auto_differentiation(c(1), f)
  expect_equal(result[1, 1], cos(1))
})

test_that("autodiff derivative of cos(x)", {
  f <- function(x) cos(x[1])
  result <- auto_differentiation(c(1), f)
  expect_equal(result[1, 1], -sin(1))
})

test_that("autodiff derivative of tan(x)", {
  f <- function(x) tan(x[1])
  result <- auto_differentiation(c(0.5), f)
  expect_equal(result[1, 1], 1 / cos(0.5)^2)
})

# ---- chain rule ----

test_that("autodiff chain rule: exp(x^2)", {
  f <- function(x) exp(x[1]^2)
  result <- auto_differentiation(c(2), f)
  # d/dx exp(x^2) = 2x * exp(x^2) = 4 * exp(4)
  expect_equal(result[1, 1], 4 * exp(4))
})

test_that("autodiff chain rule: log(1 + exp(x))", {
  f <- function(x) log(1 + exp(x[1]))
  result <- auto_differentiation(c(0), f)
  # d/dx log(1 + exp(x)) = exp(x) / (1 + exp(x)) = 0.5
  expect_equal(result[1, 1], 0.5)
})

test_that("autodiff chain rule: sin(exp(x))", {
  f <- function(x) sin(exp(x[1]))
  result <- auto_differentiation(c(0), f)
  # d/dx sin(exp(x)) = exp(x) * cos(exp(x)) = cos(1)
  expect_equal(result[1, 1], cos(1))
})

# ---- negation ----

test_that("autodiff negation", {
  f <- function(x) -x[1]^2
  result <- auto_differentiation(c(3), f)
  expect_equal(result[1, 1], -6)
})

# ---- division ----

test_that("autodiff quotient rule: x1 / x2", {
  f <- function(x) x[1] / x[2]
  result <- auto_differentiation(c(6, 3), f)
  # df/dx1 = 1/x2 = 1/3
  # df/dx2 = -x1/x2^2 = -6/9 = -2/3
  expect_equal(result[1, 1], 1 / 3)
  expect_equal(result[1, 2], -6 / 9)
})

test_that("autodiff reverse division: constant / x", {
  f <- function(x) 1 / x[1]
  result <- auto_differentiation(c(2), f)
  # d/dx (1/x) = -1/x^2 = -0.25
  expect_equal(result[1, 1], -0.25)
})

# ---- power rule ----

test_that("autodiff constant^x", {
  f <- function(x) 2^x[1]
  result <- auto_differentiation(c(3), f)
  # d/dx 2^x = 2^x * ln(2) = 8 * ln(2)
  expect_equal(result[1, 1], 8 * log(2))
})

# ---- comparisons ----

test_that("autodiff handles comparison operators", {
  # Comparisons return plain logical for if/else control flow
  f <- function(x) {
    if (x[1] > 0) x[1]^2 else -x[1]^2
  }
  result_pos <- auto_differentiation(c(3), f)
  expect_equal(result_pos[1, 1], 6) # 2x = 6

  result_neg <- auto_differentiation(c(-3), f)
  expect_equal(result_neg[1, 1], 6) # -2x = 6
})

# ---- abs ----

test_that("autodiff absolute value", {
  f <- function(x) abs(x[1])
  result_pos <- auto_differentiation(c(3), f)
  expect_equal(result_pos[1, 1], 1)

  result_neg <- auto_differentiation(c(-3), f)
  expect_equal(result_neg[1, 1], -1)
})

# ---- sum ----

test_that("autodiff derivative through sum", {
  # f(theta) = sum(y - theta) where y is data
  y <- c(1, 2, 3, 4, 5)
  f <- function(x) sum(y - x[1])
  result <- auto_differentiation(c(0), f)
  # d/dtheta sum(y - theta) = -n = -5
  expect_equal(result[1, 1], -5)
})

test_that("autodiff derivative through weighted sum", {
  w <- c(1, 2, 3)
  y <- c(10, 20, 30)
  f <- function(x) sum(w * (y - x[1]))
  result <- auto_differentiation(c(0), f)
  # d/dtheta sum(w * (y - theta)) = -sum(w) = -6
  expect_equal(result[1, 1], -6)
})

# ---- matches numerical differentiation ----

test_that("autodiff matches approx_differentiation for complex function", {
  f <- function(x) {
    c(x[1]^2 + sin(x[2]), exp(x[1]) * x[2])
  }

  theta <- c(1.5, 0.8)
  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta)

  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("autodiff matches approx for regression-like function", {
  # Summed regression EEs using element-wise operations (no matrix ops)
  set.seed(42)
  n <- 20
  x <- rnorm(n)
  y <- rnorm(n)

  f <- function(theta) {
    resid <- y - (theta[1] + theta[2] * x)
    c(sum(resid), sum(resid * x))
  }

  theta <- c(0.5, 1.2)
  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta)

  expect_equal(exact, approx, tolerance = 1e-5)
})

# ---- integration with MEstimator ----

test_that("autodiff gives same variance as approx in MEstimator", {
  set.seed(42)
  n <- 100
  y <- rnorm(n, mean = 5)

  # Custom psi using element-wise ops (no matrix() calls)
  psi <- function(theta) {
    y - theta[1]
  }

  # With numerical differentiation (default)
  m_approx <- MEstimator(stacked_equations = psi, init = 0)
  m_approx <- estimate(m_approx)

  # With exact autodiff
  m_exact <- MEstimator(stacked_equations = psi, init = 0)
  m_exact <- estimate(m_exact, deriv_method = "exact")

  expect_equal(m_exact@theta, m_approx@theta)
  expect_equal(m_exact@variance, m_approx@variance, tolerance = 1e-6)
})

# ---- matrix operations on tangent-carrying values ---------------------------
#
# Exact forward-mode autodiff must support the matrix operations that every
# built-in estimating equation relies on: `%*%`, `matrix()`, `rbind()`,
# `cbind()`, `t()`, and multi-dimensional indexing on values that carry
# tangents (the parameter vector and the results of `%*%`). Without these, no
# built-in `ee_*` function works in exact mode, which is the package-wide gap
# tracked by this epic. These matrix operations are implemented in the autodiff
# system, so the tests below pass.
#
# Each derivative is cross-checked two ways:
#   1. against `approx_differentiation()` with method "capprox" (the central
#      difference approximation) on the *same* function, and
#   2. against a hand-derived analytic Jacobian, so the assertion discriminates
#      a correct implementation from a wrong one (a comparison against capprox
#      alone can pass vacuously when the true Jacobian is degenerate).
#
# The capprox cross-check uses tolerance 1e-5 because the central difference
# approximation carries roughly 1e-6 numerical error at the default step size
# (dominated by floating-point cancellation), consistent with the existing
# autodiff-vs-approx tests above. The analytic cross-check uses 1e-8 because
# exact autodiff should reproduce the closed-form Jacobian to near machine
# precision.
#
# Each `f` below exercises a matrix operation on tangent-carrying values; the
# reductions are written so that each `f` returns a plain numeric vector under a
# numeric argument, which keeps the capprox cross-check valid.

test_that("autodiff differentiates matrix-vector product X %*% theta", {
  # Design-matrix pattern: constant X times the parameter vector.
  set.seed(2024)
  X <- matrix(rnorm(6), nrow = 3, ncol = 2) # constant design matrix (data)
  theta <- c(0.5, -1.2)

  f <- function(th) {
    eta <- X %*% th # length-3 predictor carrying tangents
    c(eta[1], eta[2], eta[3]) # reduce to a plain vector by indexing
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Jacobian of X %*% theta with respect to theta is exactly X.
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, X, tolerance = 1e-8)
})

test_that("autodiff differentiates a regression-score EE using X %*% theta", {
  # Realistic estimating-equation shape and the package-wide blocking case:
  # the summed linear-regression score t(X) %*% (y - X %*% theta), whose
  # Jacobian (the bread up to sign) is -X'X.
  set.seed(7)
  X <- matrix(rnorm(18), nrow = 6, ncol = 3) # constant design matrix (data)
  y <- rnorm(6) # constant outcomes (data)
  theta <- c(0.3, -0.8, 1.1)

  f <- function(th) {
    resid <- y - X %*% th # length-6 residuals carrying tangents
    sc <- t(X) %*% resid # length-3 score carrying tangents
    c(sc[1], sc[2], sc[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, -t(X) %*% X, tolerance = 1e-8)
})

test_that("autodiff differentiates matrix-matrix product A %*% Theta", {
  # Constant matrix times a matrix assembled from tangent-carrying parameters.
  A <- matrix(c(1, 0, 2, 3), nrow = 2) # constant, [[1, 2], [0, 3]]
  theta <- c(0.5, -1.2, 2.0, 0.3)

  f <- function(x) {
    Theta <- matrix(c(x[1], x[2], x[3], x[4]), nrow = 2)
    P <- A %*% Theta # 2-by-2 product carrying tangents
    c(P[1, 1], P[2, 1], P[1, 2], P[2, 2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1 + 2 x2, 3 x2, x3 + 2 x4, 3 x4).
  analytic <- rbind(
    c(1, 2, 0, 0),
    c(0, 3, 0, 0),
    c(0, 0, 1, 2),
    c(0, 0, 0, 3)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates vector-matrix product theta %*% M", {
  # Tangent-carrying vector times a constant matrix.
  M <- matrix(c(1, 0, 2, 0, 1, 3, 2, 1, 0), nrow = 3, byrow = TRUE) # constant
  theta <- c(0.5, -1.2, 2.0)

  f <- function(x) {
    r <- x %*% M # length-3 result carrying tangents
    c(r[1], r[2], r[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # d(sum_j x_j M[j, k]) / d x_j = M[j, k], so the Jacobian is t(M).
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, t(M), tolerance = 1e-8)
})

test_that("autodiff differentiates through matrix() construction from tangents", {
  # Reshape the tangent-carrying parameter vector into a matrix, then use it.
  theta <- c(0.7, -0.4, 1.3, 0.9)

  f <- function(x) {
    M <- matrix(x, nrow = 2) # column-major reshape of a tangent vector
    c(M[1, 1]^2, M[2, 1], M[1, 2], M[2, 2]^3)
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1^2, x2, x3, x4^3).
  analytic <- rbind(
    c(2 * theta[1], 0, 0, 0),
    c(0, 1, 0, 0),
    c(0, 0, 1, 0),
    c(0, 0, 0, 3 * theta[4]^2)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates through rbind() of tangent-carrying rows", {
  theta <- c(0.5, -1.2, 2.0)

  f <- function(x) {
    a <- x[1:2] # tangent-carrying vector
    b <- x[2:3] # tangent-carrying vector
    M <- rbind(a, b) # 2-by-2 stacked by rows
    c(M[1, 1]^2, M[1, 2], M[2, 1], M[2, 2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1^2, x2, x2, x3).
  analytic <- rbind(
    c(2 * theta[1], 0, 0),
    c(0, 1, 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates through cbind() of tangent-carrying columns", {
  theta <- c(0.5, -1.2, 2.0)

  f <- function(x) {
    a <- x[1:2] # tangent-carrying vector
    b <- x[2:3] # tangent-carrying vector
    M <- cbind(a, b) # 2-by-2 stacked by columns
    c(M[1, 1], M[1, 2]^2, M[2, 1], M[2, 2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1, x2^2, x2, x3).
  analytic <- rbind(
    c(1, 0, 0),
    c(0, 2 * theta[2], 0),
    c(0, 1, 0),
    c(0, 0, 1)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates through t() transpose of a tangent matrix", {
  theta <- c(0.7, -0.4, 1.3, 0.9)

  f <- function(x) {
    M <- matrix(x, nrow = 2) # 2-by-2 tangent-carrying matrix
    Mt <- t(M) # transpose
    c(Mt[1, 1], Mt[1, 2]^2, Mt[2, 1], Mt[2, 2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Transpose maps (M[1,1], M[2,1], M[1,2], M[2,2]) = (x1, x2, x3, x4) to
  # (Mt[1,1], Mt[1,2], Mt[2,1], Mt[2,2]) = (x1, x2, x3, x4); outputs apply the
  # square to Mt[1,2] = x2, giving (x1, x2^2, x3, x4).
  analytic <- rbind(
    c(1, 0, 0, 0),
    c(0, 2 * theta[2], 0, 0),
    c(0, 0, 1, 0),
    c(0, 0, 0, 1)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff supports [, [[, and [i, j] indexing of tangent values", {
  theta <- c(1.5, -0.5, 2.0, 0.7)

  f <- function(x) {
    v <- x[1:3] # `[` slice of a tangent-carrying vector
    M <- matrix(x, nrow = 2) # 2-by-2 tangent-carrying matrix
    c(v[[1]]^2, v[2], M[2, 1] * M[1, 2]) # `[[`, `[`, and `[i, j]` indexing
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1^2, x2, x2 * x3).
  analytic <- rbind(
    c(2 * theta[1], 0, 0, 0),
    c(0, 1, 0, 0),
    c(0, theta[3], theta[2], 0)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff selects whole rows and columns of a tangent matrix", {
  # M[i, ] and M[, j] must return the row/column, not the (i, 1) or (1, j)
  # scalar. The asymmetric 2-by-3 shape discriminates a row/column mix-up.
  theta <- c(0.7, -0.4, 1.3, 0.9, 2.1, -1.5)

  f <- function(x) {
    M <- matrix(x, nrow = 2) # 2-by-3, column-major
    row1 <- M[1, ] # whole first row -> (x1, x3, x5)
    col2 <- M[, 2] # whole second column -> (x3, x4)
    c(row1[1], row1[2]^2, row1[3], col2[1], col2[2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Outputs are (x1, x3^2, x5, x3, x4).
  analytic <- rbind(
    c(1, 0, 0, 0, 0, 0),
    c(0, 0, 2 * theta[3], 0, 0, 0),
    c(0, 0, 0, 0, 1, 0),
    c(0, 0, 1, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("M[i, ] returns the whole row rather than the leading scalar", {
  # Direct guard: on a 2-by-3 tangent matrix, M[1, ] must carry three elements.
  pta <- primal_tangent_array(matrix(1:6, nrow = 2), matrix(7:12, nrow = 2))

  row1 <- pta[1, ]
  expect_s3_class(row1, "PrimalTangentArray")
  expect_equal(length(row1), 3)
  expect_equal(as.vector(row1$primal), c(1, 3, 5))
  expect_equal(as.vector(row1$tangent), c(7, 9, 11))

  col2 <- pta[, 2]
  expect_equal(as.vector(col2$primal), c(3, 4))
  expect_equal(as.vector(col2$tangent), c(9, 10))
})

# ---- distribution and gamma-family tangent rules -----------------------------
#
# Exact forward-mode autodiff must carry tangents through the distribution and
# gamma-family functions that the built-in estimating equations rely on: the
# logistic transform (`inverse_logit`, whose kernel is `plogis`), the standard
# normal CDF and PDF (`standard_normal_cdf`/`pnorm`, `standard_normal_pdf`/
# `dnorm`), the polygamma family (`deli_polygamma`, `deli_digamma`), the base R
# Math-group members `lgamma`, `trigamma`, `log1p`, and `expm1`, and the running
# reduction `cumsum`. Python's `PrimalTangentPairs` implements each of these
# (delicatessen/derivative.py); the R side dispatches through the same surface
# functions (delicatessen/utilities.py) and through R's Math group generic.
# The wrappers dispatch on tangent-carrying inputs and the Math group generics
# carry the tangent rules, so the tests below pass.
#
# Each derivative is cross-checked two ways, matching the matrix-operation block
# above:
#   1. against `approx_differentiation()` with method "capprox" (central
#      difference) on the *same* function, and
#   2. against the hand-derived closed-form derivative, so the assertion
#      discriminates a correct rule from a wrong one rather than passing
#      vacuously against a degenerate approximation.
#
# The capprox cross-check uses tolerance 1e-5 because the central difference
# approximation carries roughly 1e-6 numerical error at the default step size
# (dominated by floating-point cancellation), consistent with the autodiff-vs-
# approx tests above. The analytic cross-check uses 1e-8 because exact autodiff
# should reproduce the closed form to near machine precision.
#
# The closed forms exercised below (verified independently against capprox):
#   d/dx plogis(x)      = plogis(x) * (1 - plogis(x))
#   d/dx pnorm(x)       = dnorm(x)
#   d/dx dnorm(x)       = -x * dnorm(x)
#   d/dx polygamma(n,x) = polygamma(n + 1, x)   (psigamma(x, deriv = n + 1))
#   d/dx digamma(x)     = trigamma(x)           (psigamma(x, deriv = 1))
#   d/dx lgamma(x)      = digamma(x)
#   d/dx trigamma(x)    = psigamma(x, deriv = 2)
#   d/dx log1p(x)       = 1 / (1 + x)
#   d/dx expm1(x)       = exp(x)
#   Jacobian of cumsum  = lower-triangular ones
#
# Each `f` returns a plain numeric vector under a numeric argument, which keeps
# the capprox cross-check valid.

# ---- scalar PrimalTangent: distribution functions ---------------------------

test_that("autodiff differentiates inverse_logit (plogis kernel)", {
  x0 <- 0.4
  f <- function(x) inverse_logit(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  analytic <- plogis(x0) * (1 - plogis(x0))
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates standard_normal_cdf (pnorm)", {
  x0 <- 0.5
  f <- function(x) standard_normal_cdf(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx pnorm(x) = dnorm(x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], dnorm(x0), tolerance = 1e-8)
})

test_that("autodiff differentiates standard_normal_pdf (dnorm)", {
  x0 <- 0.5
  f <- function(x) standard_normal_pdf(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx dnorm(x) = -x * dnorm(x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], -x0 * dnorm(x0), tolerance = 1e-8)
})

test_that("autodiff differentiates standard_normal_cdf under the chain rule", {
  # Guards that the normal-CDF rule composes: d/dx pnorm(x^2) = dnorm(x^2)*2x.
  x0 <- 0.7
  f <- function(x) standard_normal_cdf(x[1]^2)

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  analytic <- dnorm(x0^2) * 2 * x0
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
})

# ---- scalar PrimalTangent: polygamma family ---------------------------------

test_that("autodiff differentiates deli_polygamma at orders 0, 1, and 2", {
  x0 <- 2.5
  # d/dx polygamma(n, x) = polygamma(n + 1, x)
  for (n in 0:2) {
    f <- function(x) deli_polygamma(n, x[1])

    exact <- auto_differentiation(x0, f)
    approx <- approx_differentiation(f, x0, method = "capprox")

    analytic <- psigamma(x0, deriv = n + 1)
    expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
    expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
  }
})

test_that("autodiff differentiates deli_digamma", {
  x0 <- 2.5
  f <- function(x) deli_digamma(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx digamma(x) = trigamma(x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], trigamma(x0), tolerance = 1e-8)
})

# ---- scalar PrimalTangent: base R Math-group members ------------------------

test_that("autodiff differentiates lgamma", {
  x0 <- 2.5
  f <- function(x) lgamma(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx lgamma(x) = digamma(x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], digamma(x0), tolerance = 1e-8)
})

test_that("autodiff differentiates trigamma", {
  x0 <- 2.5
  f <- function(x) trigamma(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx trigamma(x) = psigamma(x, deriv = 2)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], psigamma(x0, deriv = 2), tolerance = 1e-8)
})

test_that("autodiff differentiates log1p", {
  x0 <- 0.5
  f <- function(x) log1p(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx log1p(x) = 1 / (1 + x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], 1 / (1 + x0), tolerance = 1e-8)
})

test_that("autodiff differentiates expm1", {
  x0 <- 0.5
  f <- function(x) expm1(x[1])

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  # d/dx expm1(x) = exp(x)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], exp(x0), tolerance = 1e-8)
})

# ---- PrimalTangentVector: partials route to the correct columns -------------

test_that("autodiff routes gamma-family partials across a PrimalTangentVector", {
  # Two independent parameters through two different rules; the Jacobian must be
  # diagonal, which discriminates a column mix-up from a correct routing.
  theta <- c(2.5, 3.1)

  f <- function(x) c(deli_digamma(x[1]), lgamma(x[2]))

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(
    c(trigamma(theta[1]), 0),
    c(0, digamma(theta[2]))
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff routes distribution partials across a PrimalTangentVector", {
  theta <- c(0.3, -0.6)

  f <- function(x) c(standard_normal_cdf(x[1]), inverse_logit(x[2]))

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(
    c(dnorm(theta[1]), 0),
    c(0, plogis(theta[2]) * (1 - plogis(theta[2])))
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

# ---- tangent vectors (PrimalTangentArray): cumsum ---------------------------

test_that("autodiff differentiates cumsum over a tangent-carrying vector", {
  # cumsum is linear, so its Jacobian is the lower-triangular matrix of ones.
  theta <- c(1.0, 2.0, 3.0)

  f <- function(x) {
    v <- matrix(x, ncol = 1) # tangent-carrying column vector
    s <- cumsum(v) # running total, one output per element
    c(s[1], s[2], s[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(
    c(1, 0, 0),
    c(1, 1, 0),
    c(1, 1, 1)
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates cumsum composed with a nonlinear map", {
  # cumsum(x^2): output k is sum_{j<=k} x_j^2, so the Jacobian is lower
  # triangular with entry (k, j) = 2 x_j for j <= k. The nonlinear entries
  # discriminate the running-sum structure from a plain identity.
  theta <- c(0.7, -1.1, 1.4)

  f <- function(x) {
    v <- matrix(x, ncol = 1)
    s <- cumsum(v^2)
    c(s[1], s[2], s[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  g <- 2 * theta
  analytic <- rbind(
    c(g[1], 0, 0),
    c(g[1], g[2], 0),
    c(g[1], g[2], g[3])
  )
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

# ---- unwired array Math members abort with a rendered message ----------------

test_that("a genuinely unwired Math member aborts on a tangent array with a rendered message", {
  # `gamma` is not part of the scalar tangent surface, so it must reach the abort
  # default on a tangent array too. This pins the cli message rendering:
  # interpolating `.Generic` directly would raise "Invalid cli literal:
  # `{.Generic}` starts with a dot" instead of the designed message.
  pta <- primal_tangent_array(matrix(1:4, nrow = 2), matrix(1:4, nrow = 2))
  expect_error(
    gamma(pta),
    "Math function.*gamma.*not supported for tangent arrays"
  )
})

# ---- EE-shaped composite: logistic score with plogis(X %*% theta) -----------

test_that("autodiff differentiates a logistic score using plogis(X %*% theta)", {
  # Realistic estimating-equation shape: the summed logistic-regression score
  # t(X) %*% (y - plogis(X %*% theta)), which now routes a tangent-carrying
  # predictor through inverse_logit (the plogis kernel) elementwise. Its
  # Jacobian is the negative weighted cross-product -X' W X with
  # W = diag(p (1 - p)) and p = plogis(X %*% theta).
  set.seed(11)
  n <- 8
  X <- cbind(1, rnorm(n), rnorm(n)) # constant design matrix (data)
  y <- rbinom(n, 1, 0.5) # constant outcomes (data)
  theta <- c(0.2, -0.5, 0.3)

  f <- function(th) {
    eta <- X %*% th # length-n predictor carrying tangents
    p <- inverse_logit(eta) # elementwise logistic transform
    sc <- t(X) %*% (y - p) # length-3 score carrying tangents
    c(sc[1], sc[2], sc[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  pv <- as.numeric(plogis(X %*% theta))
  analytic <- -t(X) %*% (X * (pv * (1 - pv)))
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

# ---- vector-primal Math operations and prod ---------------------------------
#
# Exact forward-mode autodiff must carry tangents through `abs`, `floor`, and
# `ceiling` when the primal is a vector rather than a scalar, and must support
# `prod` as a differentiable reduction. The vector-primal shape arises whenever a
# tangent-carrying scalar broadcasts over a data vector (for example
# `abs(x[1] - y)` yields a PrimalTangent whose primal is a vector and whose
# tangent is the broadcast scalar) and whenever a tangent-carrying vector or
# matrix flows through `%*%`, `matrix()`, or indexing (a PrimalTangentArray).
#
# The current `abs`/`floor`/`ceiling` branches select with a scalar `if (p > 0)`
# style test, which under R >= 4.2 raises "the condition has length > 1" on a
# vector primal (confirmed as the observed failure mode on R 4.6). `prod` is
# unsupported and stops explicitly for a PrimalTangent, and reaches an
# unclassed-list error for a PrimalTangentArray. Every test below therefore
# fails today, either by the length > 1 error, the explicit `prod` stop, or the
# Math-array abort, and pins the intended post-fix behavior.
#
# Conventions mirrored from Python's PrimalTangentPairs (delicatessen/
# derivative.py):
#   __abs__  : primal = |p|; tangent = sign(p) * t elementwise, i.e. +t where
#              p > 0 and -t where p < 0. At exactly p == 0 the derivative is
#              undefined and the tangent is NaN (the scalar branch already
#              returns NaN there; the vector rule preserves that).
#   __floor__: primal = floor(p); tangent = 0 away from integers (p %% 1 != 0)
#              and NaN at integers (p %% 1 == 0), where the derivative is
#              undefined. This is the floor/ceiling convention adopted here.
#   __ceil__ : primal = ceiling(p); tangent = 0 away from integers, NaN at
#              integers, identical convention to floor.
#   prod     : d prod(x) / dx_i = prod(x[-i]) via the product rule. A built-in
#              estimating equation search found no use of prod(); a scalar/vector
#              test therefore suffices and no PrimalTangentArray prod is wired.
#
# `prod(c(x[1], x[2], ...))` cannot be used because `c.PrimalTangent` returns an
# unclassed list, so the reduction does not dispatch to the Summary generic (the
# same limitation the existing `sum` tests avoid by summing a single vector-
# primal pair). The natural multiple-argument form `prod(x[1], x[2], ...)`
# dispatches to Summary.PrimalTangent and is used instead.
#
# Cross-checks follow the matrix and gamma-family blocks above: capprox at tolerance 1e-5
# (the central difference carries roughly 1e-6 numerical error at the default
# step size) and the hand-derived analytic derivative at 1e-8 (exact autodiff
# reproduces the closed form to near machine precision). Evaluation points for
# the capprox checks sit away from the kink at 0 and away from integer
# boundaries so the central difference is valid; the direct guards pin the
# behavior exactly at those boundaries.

# ---- abs on vector primals --------------------------------------------------

test_that("autodiff differentiates abs over a vector primal (PrimalTangent)", {
  # `x[1] - y` is a PrimalTangent with a vector primal and a broadcast scalar
  # tangent; abs must apply sign(p) * t elementwise so the summed derivative is
  # sum(sign(x1 - y)). The primal entries sit away from 0 for the capprox check.
  y <- c(-2, -1, 3)
  x0 <- 0.5 # x1 - y = c(2.5, 1.5, -2.5): mixed sign
  f <- function(x) sum(abs(x[1] - y))

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  analytic <- sum(sign(x0 - y)) # 1 + 1 - 1 = 1
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
})

test_that("abs on a vector-primal PrimalTangent applies sign(p) * t elementwise", {
  # Direct guard covering the mixed-sign case the scalar if() cannot handle,
  # including the p == 0 kink where the tangent is NaN (Python convention).
  pt <- primal_tangent(c(-2, 0, 3), 1) # vector primal, scalar broadcast tangent
  r <- abs(pt)

  expect_equal(r$primal, c(2, 0, 3))
  expect_equal(r$tangent[1], -1) # sign(-2) * 1
  expect_true(is.nan(r$tangent[2])) # undefined at 0
  expect_equal(r$tangent[3], 1) # sign(3) * 1
})

test_that("autodiff differentiates abs over a PrimalTangentArray", {
  # `X %*% theta` is a PrimalTangentArray; abs must apply elementwise so the
  # Jacobian of abs(X theta) with respect to theta is diag(sign(eta)) %*% X.
  X <- matrix(c(1, -1, 0.5, 2, 1, -1), nrow = 3, ncol = 2) # constant (data)
  theta <- c(0.9, -0.8) # eta = c(-0.7, -1.7, 1.25): away from 0

  f <- function(th) {
    eta <- X %*% th # length-3 predictor carrying tangents
    a <- abs(eta) # elementwise abs on the tangent array
    c(a[1], a[2], a[3]) # reduce to a plain vector by indexing
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- sign(eta) * X
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("abs on a PrimalTangentArray applies sign(p) * t elementwise", {
  # Direct guard on a mixed-sign matrix; both slots keep their 2-by-2 shape.
  pta <- primal_tangent_array(
    matrix(c(-1, 2, -3, 4), nrow = 2),
    matrix(c(1, 1, 1, 1), nrow = 2)
  )
  r <- abs(pta)

  expect_s3_class(r, "PrimalTangentArray")
  expect_equal(r$primal, matrix(c(1, 2, 3, 4), nrow = 2))
  expect_equal(r$tangent, matrix(c(-1, 1, -1, 1), nrow = 2)) # sign(p) * t
})

test_that("abs on a PrimalTangentVector applies sign(p) * t per element", {
  # A slice or the whole parameter vector is a PrimalTangentVector; abs must
  # preserve each element's independent tangent. Distinct tangents discriminate
  # a correct per-element rule from a scalar-broadcast collapse. `pt_arrays`
  # normalizes whatever container the fix returns into parallel numeric vectors.
  ptv <- primal_tangent_vector(list(
    primal_tangent(-2, 1),
    primal_tangent(3, 2),
    primal_tangent(-0.5, 3)
  ))
  parts <- pt_arrays(abs(ptv))

  expect_equal(parts$primal, c(2, 3, 0.5))
  expect_equal(parts$tangent, c(-1, 2, -3)) # sign(p) * t, per element
})

# ---- floor and ceiling on vector primals ------------------------------------

test_that("autodiff differentiates floor over a vector primal (zero off integers)", {
  # floor has zero derivative away from integers; the Jacobian is 0. The value
  # of the test is confirming the vector primal flows through without the
  # length > 1 error and returns a zero (not NaN) tangent off the integers.
  y <- c(0.3, 1.7, 2.2)
  x0 <- 0.5 # x1 + y = c(0.8, 2.2, 2.7): non-integer
  f <- function(x) sum(floor(x[1] + y))

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], 0, tolerance = 1e-8)
})

test_that("floor on a vector-primal PrimalTangent: zero off integers, NaN at integers", {
  # Mixed integer / non-integer vector: the branch must be selected elementwise,
  # which a scalar if() cannot do. Pins the floor convention exactly.
  pt <- primal_tangent(c(2, 2.5), 1)
  r <- floor(pt)

  expect_equal(r$primal, c(2, 2))
  expect_true(is.nan(r$tangent[1])) # integer: derivative undefined
  expect_equal(r$tangent[2], 0) # non-integer: derivative 0
})

test_that("autodiff differentiates ceiling over a vector primal (zero off integers)", {
  y <- c(0.3, 1.7, 2.2)
  x0 <- 0.5 # x1 + y = c(0.8, 2.2, 2.7): non-integer
  f <- function(x) sum(ceiling(x[1] + y))

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
  expect_equal(exact[1, 1], 0, tolerance = 1e-8)
})

test_that("ceiling on a vector-primal PrimalTangent: zero off integers, NaN at integers", {
  pt <- primal_tangent(c(3, 2.5), 1)
  r <- ceiling(pt)

  expect_equal(r$primal, c(3, 3))
  expect_true(is.nan(r$tangent[1])) # integer: derivative undefined
  expect_equal(r$tangent[2], 0) # non-integer: derivative 0
})

test_that("floor and ceiling on a PrimalTangentArray follow the integer convention", {
  # Elementwise branch selection on a matrix: one integer entry (NaN tangent)
  # and three non-integer entries (zero tangent).
  pta <- primal_tangent_array(
    matrix(c(0.8, 2, -1.5, 2.7), nrow = 2), # entry [2, 1] == 2 is the integer
    matrix(c(1, 1, 1, 1), nrow = 2)
  )

  rf <- floor(pta)
  expect_s3_class(rf, "PrimalTangentArray")
  expect_equal(rf$primal, matrix(c(0, 2, -2, 2), nrow = 2))
  expect_equal(rf$tangent[1, 1], 0)
  expect_true(is.nan(rf$tangent[2, 1])) # integer entry
  expect_equal(rf$tangent[1, 2], 0)
  expect_equal(rf$tangent[2, 2], 0)

  rc <- ceiling(pta)
  expect_equal(rc$primal, matrix(c(1, 2, -1, 3), nrow = 2))
  expect_equal(rc$tangent[1, 1], 0)
  expect_true(is.nan(rc$tangent[2, 1])) # integer entry
  expect_equal(rc$tangent[1, 2], 0)
  expect_equal(rc$tangent[2, 2], 0)
})

# ---- floor and ceiling on non-finite primals --------------------------------
#
# Only a finite integer primal is a kink where the derivative is undefined and
# the tangent is NaN. Python's __floor__ and __ceil__ (delicatessen/
# derivative.py) test `self.primal % 1 == 0`, which is False for both a NaN
# primal (`nan % 1` is nan) and an infinite primal (`inf % 1` is nan), so both
# take the else branch and return a 0 tangent alongside the floor/ceiling primal
# (NaN stays NaN, Inf stays Inf). These guards pin that parity: R's
# `NaN %% 1 == 0` is NA rather than FALSE, so an unguarded rule would leak an NA
# tangent instead of the 0 Python returns.

test_that("floor and ceiling on a NaN primal give a 0 tangent (Python parity)", {
  pt <- primal_tangent(NaN, 1)

  rf <- floor(pt)
  expect_true(is.nan(rf$primal))
  expect_equal(rf$tangent, 0)

  rc <- ceiling(pt)
  expect_true(is.nan(rc$primal))
  expect_equal(rc$tangent, 0)
})

test_that("floor and ceiling on an infinite primal give a 0 tangent (Python parity)", {
  for (prim in c(Inf, -Inf)) {
    pt <- primal_tangent(prim, 1)

    rf <- floor(pt)
    expect_equal(rf$primal, prim)
    expect_equal(rf$tangent, 0)

    rc <- ceiling(pt)
    expect_equal(rc$primal, prim)
    expect_equal(rc$tangent, 0)
  }
})

test_that("floor and ceiling keep the finite-integer kink (NaN tangent)", {
  # The fix must not disturb the kink: a finite integer primal still gives a NaN
  # tangent, and a finite non-integer primal still gives a 0 tangent.
  pt_int <- primal_tangent(3, 1)
  expect_true(is.nan(floor(pt_int)$tangent))
  expect_true(is.nan(ceiling(pt_int)$tangent))

  pt_frac <- primal_tangent(2.5, 1)
  expect_equal(floor(pt_frac)$tangent, 0)
  expect_equal(ceiling(pt_frac)$tangent, 0)
})

test_that("floor and ceiling handle a mixed finite/NaN vector primal elementwise", {
  # A NaN entry alongside a finite integer and a finite non-integer: NaN and the
  # non-integer both give a 0 tangent, and only the finite integer gives NaN.
  pt <- primal_tangent(c(NaN, 2, 2.5), 1)

  rf <- floor(pt)
  expect_true(is.nan(rf$primal[1]))
  expect_equal(rf$tangent[1], 0) # NaN primal: 0 tangent
  expect_true(is.nan(rf$tangent[2])) # finite integer: kink
  expect_equal(rf$tangent[3], 0) # finite non-integer: 0 tangent
})

# ---- prod as a differentiable reduction -------------------------------------

test_that("autodiff differentiates prod of three parameters", {
  # d prod(x) / dx_i = prod(x[-i]); distinct magnitudes including a negative so
  # the gradient entries are distinct and well above the tolerance floor.
  theta <- c(2, -3, 4)
  f <- function(x) prod(x[1], x[2], x[3])

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(c(-3 * 4, 2 * 4, 2 * -3)) # c(-12, 8, -6)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates prod of four parameters", {
  theta <- c(1.5, -2, 3, 0.5)
  f <- function(x) prod(x[1], x[2], x[3], x[4])

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(vapply(
    seq_along(theta),
    function(i) prod(theta[-i]),
    numeric(1)
  )) # c(-3, 2.25, -1.5, -9)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("log aborts when the base carries a tangent under exact autodiff", {
  # log(x, base) differentiates with respect to x only; a tangent-carrying base
  # would need an extra term. Rather than silently return the constant-base
  # derivative, the exact pass must abort.
  f <- function(theta) log(theta[1], base = theta[2])
  expect_error(
    auto_differentiation(c(5, 2), f),
    "base"
  )
})

# ---- unsupported-operator messages on tangent arrays -------------------------
#
# Ops.PrimalTangentArray previously interpolated `{.val {.Generic}}` directly.
# cli rejects an interpolated expression that begins with a dot, so a user
# hitting an unwired operator saw "Invalid cli literal" instead of the intended
# message. Binding `.Generic` to a local first renders the operator correctly on
# both the binary and the unary abort paths.

test_that("an unsupported binary operator on a tangent array names the operator", {
  pta <- primal_tangent_array(
    matrix(c(1, 2, 3, 4), nrow = 2),
    matrix(c(1, 1, 1, 1), nrow = 2)
  )
  expect_error(
    pta %% pta,
    'Operator "%%" is not supported for tangent arrays.',
    fixed = TRUE
  )
})

test_that("an unsupported unary operator on a tangent array names the operator", {
  pta <- primal_tangent_array(
    matrix(c(1, 2, 3, 4), nrow = 2),
    matrix(c(1, 1, 1, 1), nrow = 2)
  )
  expect_error(
    !pta,
    'Unary "!" is not supported for tangent arrays.',
    fixed = TRUE
  )
})

# ---- elementwise Math on tangent arrays --------------------------------------
#
# Math.PrimalTangentArray previously wired only `exp`, `cumsum`, `abs`, `floor`,
# and `ceiling`; every other elementwise math member (`log`, `sqrt`, the
# trigonometric functions, the gamma family) aborted on a tangent-carrying
# vector or matrix even though the scalar surface supports them. The built-in
# estimating equations reach these members on a tangent-carrying predictor (for
# example `log(X %*% theta)`), so the array surface must mirror the scalar
# surface exactly. Each rule below is applied to the tangent-carrying predictor
# `X %*% theta`, which is a PrimalTangentArray.
#
# Cross-checks follow the established autodiff convention: capprox at tolerance
# 1e-5 (the central difference carries roughly 1e-6 numerical error at the
# default step size) and the hand-derived analytic Jacobian at 1e-8 (exact
# autodiff reproduces the closed form to near machine precision). Evaluation
# points keep the predictor inside each rule's domain (positive for `log`,
# `sqrt`, and `lgamma`). Each `f` returns a plain numeric vector under a numeric
# argument, which keeps the capprox cross-check valid.

test_that("autodiff differentiates log over a tangent array X %*% theta", {
  # d log(eta_i) / d theta_j = X[i, j] / eta_i, so the Jacobian is X scaled row
  # i by 1 / eta_i (column-major recycling divides each column by eta).
  X <- cbind(1, c(0.2, 0.4, 0.6)) # constant design matrix (data)
  theta <- c(1.5, 0.8) # eta = c(1.66, 1.82, 1.98): positive

  f <- function(th) {
    eta <- X %*% th
    l <- log(eta)
    c(l[1], l[2], l[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- X / eta
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates sqrt over a tangent array X %*% theta", {
  # d sqrt(eta_i) / d theta_j = X[i, j] / (2 sqrt(eta_i)).
  X <- cbind(1, c(0.2, 0.4, 0.6))
  theta <- c(1.5, 0.8) # eta positive

  f <- function(th) {
    eta <- X %*% th
    s <- sqrt(eta)
    c(s[1], s[2], s[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- X / (2 * sqrt(eta))
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates sin over a tangent array X %*% theta", {
  # d sin(eta_i) / d theta_j = cos(eta_i) * X[i, j].
  X <- cbind(1, c(0.3, -0.7, 1.1))
  theta <- c(0.4, 0.9)

  f <- function(th) {
    eta <- X %*% th
    s <- sin(eta)
    c(s[1], s[2], s[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- cos(eta) * X
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates lgamma over a tangent array X %*% theta", {
  # d lgamma(eta_i) / d theta_j = digamma(eta_i) * X[i, j].
  X <- cbind(1, c(0.2, 0.4, 0.6))
  theta <- c(1.5, 0.8) # eta positive, away from gamma poles

  f <- function(th) {
    eta <- X %*% th
    g <- lgamma(eta)
    c(g[1], g[2], g[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- digamma(eta) * X
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

# ---- log with an explicit base -----------------------------------
#
# `log(x, base)` must honor the base under exact autodiff for both the primal and
# the tangent: the natural-log tangent t / p is scaled by 1 / log(base). Under
# the bug the base was dropped and log10(x) differentiated as the natural log
# (0.5 instead of 1 / (2 log(10)) = 0.2171472 at x = 2). All three tangent
# surfaces are exercised (scalar pair, whole vector, tangent array).

test_that("autodiff differentiates log(x, base) honoring the base", {
  x0 <- 2
  f <- function(x) log(x[1], base = 10)

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  analytic <- 1 / (x0 * log(10)) # 0.2171472
  expect_equal(exact[1, 1], analytic, tolerance = 1e-8)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
})

test_that("autodiff forwards the base through the whole-vector log method", {
  # The exact repro: log(x, base = 10) on a PrimalTangentVector must forward the
  # base through Math.PrimalTangentVector rather than dropping it.
  x0 <- 2
  f <- function(x) log(x, base = 10)

  exact <- auto_differentiation(x0, f)
  approx <- approx_differentiation(f, x0, method = "capprox")

  expect_equal(exact[1, 1], 1 / (x0 * log(10)), tolerance = 1e-8)
  expect_equal(exact[1, 1], approx[1, 1], tolerance = 1e-5)
})

test_that("autodiff differentiates log(X %*% theta, base) honoring the base", {
  X <- cbind(1, c(0.2, 0.4, 0.6)) # constant design matrix (data)
  theta <- c(1.5, 0.8) # eta positive

  f <- function(th) {
    eta <- X %*% th
    l <- log(eta, base = 10)
    c(l[1], l[2], l[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta)
  analytic <- X / (eta * log(10))
  expect_equal(exact, analytic, tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("log with a base is a direct guard on scalar and array pairs", {
  pt <- primal_tangent(4, 1)
  r <- log(pt, base = 2)
  expect_s3_class(r, "PrimalTangent")
  expect_equal(r$primal, log(4, 2))
  expect_equal(r$tangent, 1 / (4 * log(2)))

  pta <- primal_tangent_array(c(2, 8), c(1, 1))
  ra <- log(pta, 10) # positional base
  expect_s3_class(ra, "PrimalTangentArray")
  expect_equal(ra$primal, log(c(2, 8), 10))
  expect_equal(ra$tangent, c(1, 1) / (c(2, 8) * log(10)))
})

test_that("log without a base still differentiates as the natural log", {
  pt <- primal_tangent(3, 1)
  r <- log(pt)
  expect_equal(r$primal, log(3))
  expect_equal(r$tangent, 1 / 3)
})

test_that("Math on a tangent array reuses the scalar rule (log direct guard)", {
  # Direct guard: log on a 2-by-2 tangent matrix applies t / p elementwise and
  # keeps the array shape, matching the scalar rule element for element.
  p <- matrix(c(1, 2, 4, 8), nrow = 2)
  t <- matrix(c(1, 1, 1, 1), nrow = 2)
  pta <- primal_tangent_array(p, t)

  r <- log(pta)
  expect_s3_class(r, "PrimalTangentArray")
  expect_equal(r$primal, log(p))
  expect_equal(r$tangent, t / p)
})

test_that("autodiff differentiates a Poisson-style score using exp(X %*% theta)", {
  # EE-shaped composite: the summed Poisson score t(X) %*% (y - exp(X %*% theta)),
  # whose Jacobian is the negative weighted cross-product -X' diag(exp(eta)) X.
  set.seed(13)
  n <- 8
  X <- cbind(1, rnorm(n), rnorm(n)) # constant design matrix (data)
  y <- rpois(n, 2) # constant outcomes (data)
  theta <- c(0.1, -0.3, 0.2)

  f <- function(th) {
    eta <- X %*% th
    mu <- exp(eta)
    sc <- t(X) %*% (y - mu)
    c(sc[1], sc[2], sc[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  mu <- as.numeric(exp(X %*% theta))
  analytic <- -t(X) %*% (X * mu)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

# ---- mixed PrimalTangentArray / scalar PrimalTangent Ops ---------------------
#
# Mixing a tangent-carrying vector or matrix with a scalar tangent pair (for
# example `(X %*% beta) / theta[k]`, the tobit and AFT scale idiom) previously
# triggered R's "Incompatible methods" fallback: the array operand selected
# `Ops.PrimalTangentArray` while the scalar operand selected `Ops.PrimalTangent`,
# and because the two were distinct functions R refused to dispatch and errored.
# The two methods are now the same function object, so either operand routes the
# call to the shared implementation. Both operand orders are exercised.
#
# The composite tests use the closed-form Jacobian at 1e-8 and capprox at 1e-5;
# the direct guards pin the primal and tangent slots exactly.

test_that("autodiff differentiates (X %*% beta) / theta[k] (tobit/AFT scale idiom)", {
  # A tangent-carrying predictor divided by a scalar tangent parameter. The
  # first two parameters form beta; the third is the scale in the denominator.
  set.seed(21)
  X <- matrix(rnorm(6), nrow = 3, ncol = 2) # constant design matrix (data)
  theta <- c(0.5, -1.2, 2.0) # beta = theta[1:2], scale = theta[3]

  f <- function(th) {
    eta <- X %*% th[1:2] # PrimalTangentArray predictor
    s <- th[3] # scalar PrimalTangent denominator
    z <- eta / s # mixed array / scalar dispatch
    c(z[1], z[2], z[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta[1:2])
  s <- theta[3]
  # d z_i / d beta_j = X[i, j] / s; d z_i / d s = -eta_i / s^2.
  analytic <- cbind(X / s, -eta / s^2)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("autodiff differentiates theta[k] * (X %*% beta) (scalar times array)", {
  set.seed(22)
  X <- matrix(rnorm(6), nrow = 3, ncol = 2) # constant design matrix (data)
  theta <- c(0.5, -1.2, 2.0) # beta = theta[1:2], scale = theta[3]

  f <- function(th) {
    eta <- X %*% th[1:2]
    s <- th[3]
    z <- s * eta # scalar * array dispatch
    c(z[1], z[2], z[3])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  eta <- as.numeric(X %*% theta[1:2])
  s <- theta[3]
  # d z_i / d beta_j = s * X[i, j]; d z_i / d s = eta_i.
  analytic <- cbind(s * X, eta)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("mixed array / scalar Ops route to the array method in both orders", {
  # Direct guards: a length-2 tangent vector meeting a scalar tangent pair. Each
  # operator and both operand orders must return a PrimalTangentArray with the
  # elementwise-broadcast primal and the correct tangent rule.
  pta <- primal_tangent_array(c(2, 4), c(1, 0)) # vector primal, tangent
  pts <- primal_tangent(2, 0.5) # scalar primal, tangent

  # division: array / scalar
  r <- pta / pts
  expect_s3_class(r, "PrimalTangentArray")
  expect_equal(r$primal, c(1, 2))
  # (t1 * p2 - p1 * t2) / p2^2 = (c(1,0)*2 - c(2,4)*0.5) / 4 = c(0.25, -0.5)
  expect_equal(r$tangent, c(0.25, -0.5))

  # division: scalar / array
  r <- pts / pta
  expect_s3_class(r, "PrimalTangentArray")
  expect_equal(r$primal, c(1, 0.5))
  # (t1 * p2 - p1 * t2) / p2^2 = (0.5*c(2,4) - 2*c(1,0)) / c(4,16)
  expect_equal(r$tangent, (0.5 * c(2, 4) - 2 * c(1, 0)) / c(4, 16))

  # multiplication both orders
  r <- pta * pts
  expect_equal(r$primal, c(4, 8))
  # t1 * p2 + p1 * t2 = c(1,0)*2 + c(2,4)*0.5 = c(3, 2)
  expect_equal(r$tangent, c(3, 2))
  r2 <- pts * pta
  expect_equal(r2$primal, c(4, 8))
  expect_equal(r2$tangent, c(3, 2))

  # subtraction both orders
  r <- pta - pts
  expect_equal(r$primal, c(0, 2))
  expect_equal(r$tangent, c(0.5, -0.5)) # t1 - t2 = c(1,0) - 0.5
  r <- pts - pta
  expect_equal(r$primal, c(0, -2))
  expect_equal(r$tangent, c(-0.5, 0.5)) # t2 - t1 = 0.5 - c(1,0)

  # addition
  r <- pta + pts
  expect_equal(r$primal, c(4, 6))
  expect_equal(r$tangent, c(1.5, 0.5)) # t1 + t2 = c(1,0) + 0.5
})

test_that("scalar / scalar PrimalTangent Ops still return a scalar pair", {
  # The shared Ops implementation must not promote a pure scalar-scalar operation
  # to an array: two scalar pairs still yield a PrimalTangent, not a
  # PrimalTangentArray.
  a <- primal_tangent(3, 1)
  b <- primal_tangent(2, 0)
  r <- a / b
  expect_s3_class(r, "PrimalTangent")
  expect_false(is_pt_array(r))
  expect_equal(r$primal, 1.5)
  expect_equal(r$tangent, (1 * 2 - 3 * 0) / 2^2)
})

# ---- pt_flatten normalization -----------------------------------------------
#
# A masked matrix()/rbind()/cbind() applied to a single scalar pair whose primal
# is a length-n vector must flatten the primal and tangent slots to matching
# lengths. This happens for a built-in EE such as ee_mean, whose
# matrix(w * (y - theta[1]), nrow = 1) hands matrix() a bare PrimalTangent with a
# length-n vector primal. When the tangent is a scalar broadcast, it must be
# recycled up to the primal length before flattening.

test_that("pt_flatten flattens a scalar pair with a vector primal", {
  # tangent already the same length as the primal (the ee_mean shape, where
  # y - theta[1] carries a length-n zero data tangent minus a scalar tangent)
  pt <- primal_tangent(c(10, 20, 30), c(-1, -1, -1))
  flat <- pt_flatten(pt)
  expect_equal(flat$primal, c(10, 20, 30))
  expect_equal(flat$tangent, c(-1, -1, -1))
})

test_that("pt_flatten recycles a scalar broadcast tangent to the primal length", {
  # a length-n primal carrying a single scalar tangent must expand element for
  # element so the tangent slot reshapes identically to the primal slot
  pt <- primal_tangent(c(2, 4, 6, 8), 0.5)
  flat <- pt_flatten(pt)
  expect_equal(flat$primal, c(2, 4, 6, 8))
  expect_equal(flat$tangent, c(0.5, 0.5, 0.5, 0.5))
})

test_that("masked matrix() reshapes a scalar pair with a vector primal", {
  # matrix() on ee_mean's PrimalTangent result
  pt <- primal_tangent(c(1, 2, 3, 4, 5), 0.5)
  m <- matrix(pt, nrow = 1)
  expect_s3_class(m, "PrimalTangentArray")
  expect_equal(m$primal, matrix(c(1, 2, 3, 4, 5), nrow = 1))
  expect_equal(m$tangent, matrix(rep(0.5, 5), nrow = 1))
})

test_that("masked rbind() binds scalar pairs with vector primals and scalar tangents", {
  # each row is a bare PrimalTangent with a length-n vector primal; a scalar
  # broadcast tangent on either row must recycle so the tangent slots bind to
  # the same 2-by-n layout as the primal slots (the ee_mean_variance shape)
  row1 <- primal_tangent(c(1, 2, 3), c(-1, -1, -1))
  row2 <- primal_tangent(c(4, 5, 6), 2)
  m <- rbind(row1, row2)
  expect_s3_class(m, "PrimalTangentArray")
  expect_equal(dim(m$primal), c(2, 3))
  expect_equal(m$tangent, rbind(c(-1, -1, -1), c(2, 2, 2)))
})

# ---- indexing a scalar PrimalTangent ----------------------------------------
#
# A scalar PrimalTangent is a single tangent-carrying value, so it indexes like
# a length-1 vector: the sole subscript selects the value itself and preserves
# the class, and any other subscript is out of bounds. The stacked delta-
# effective-dose equations (ee_emax_ed, ee_loglogistic_ed) index their scalar
# theta argument as theta[1], which under exact mode is a scalar PrimalTangent.

test_that("[ on a scalar PrimalTangent returns the pair itself for index 1", {
  x <- primal_tangent(7, 3)
  y <- x[1]
  expect_s3_class(y, "PrimalTangent")
  expect_false(is_pt_array(y))
  expect_equal(y$primal, 7)
  expect_equal(y$tangent, 3)
})

test_that("[[ on a scalar PrimalTangent returns the pair itself for index 1", {
  x <- primal_tangent(7, 3)
  y <- x[[1]]
  expect_s3_class(y, "PrimalTangent")
  expect_equal(y$primal, 7)
  expect_equal(y$tangent, 3)
})

test_that("indexing a scalar PrimalTangent out of bounds errors", {
  x <- primal_tangent(7, 3)
  expect_error(x[2], "subscript out of bounds")
  expect_error(x[[2]], "subscript out of bounds")
  expect_error(x[0], "subscript out of bounds")
})

test_that("arithmetic on an indexed scalar PrimalTangent preserves the tangent", {
  # The ee_emax_ed shape: delta/(1 - delta) * ed50 - theta[1]. With theta a
  # scalar PrimalTangent carrying a unit tangent, the result must carry -1.
  theta <- primal_tangent(40, 1)
  ed50 <- primal_tangent(10, 0)
  out <- 0.8 / (1 - 0.8) * ed50 - theta[1]
  expect_s3_class(out, "PrimalTangent")
  expect_equal(out$primal, 0.8 / 0.2 * 10 - 40)
  expect_equal(out$tangent, -1)
})

# ---- non-scalar Ops results index and transpose -------------------
#
# A scalar pair times a data vector or matrix (theta[k] * X, the interaction and
# scaling idiom) produces a PrimalTangent whose primal and tangent slots are
# themselves a vector or matrix. Selecting an element with `[`, `[[`, or `[i, j]`
# and transposing with `t()` must operate on those payloads rather than
# self-selecting the whole pair (single index) or falling through to `t.default`
# (which strips the tangent and returns NULL). The genuinely scalar case still
# self-selects for index 1.
#
# Cross-checks follow the established autodiff convention: capprox at tolerance
# 1e-5 and the hand-derived analytic Jacobian at 1e-8.

test_that("indexing a vector-primal Ops result selects a single element under exact autodiff", {
  # theta[1] * c(1, 2, 3) is a scalar PrimalTangent carrying a length-3 vector
  # primal and tangent; [1] must select the first element, not self-select the
  # whole pair (which mis-shapes the Jacobian to 3-by-3).
  theta <- c(0.1)
  Sigma <- matrix(0.01, nrow = 1, ncol = 1)
  transform <- function(th) (th[1] * c(1, 2, 3))[1]

  exact <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "exact"
  )
  capprox <- delta_method(
    theta,
    transform = transform,
    covariance = Sigma,
    deriv_method = "capprox"
  )

  expect_equal(dim(exact), c(1L, 1L))
  # d/dtheta (theta * 1) = 1, so G = [1] and G Sigma G' = 0.01.
  expect_equal(exact, matrix(0.01, nrow = 1, ncol = 1), tolerance = 1e-8)
  expect_equal(exact, capprox, tolerance = 1e-5)
})

test_that("t() differentiates a matrix-payload Ops result under exact autodiff", {
  # theta[1] * X is a scalar PrimalTangent whose slots are matrices; t() must
  # transpose both slots rather than falling through to t.default.
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2) # 2-by-3 constant (data)
  theta <- c(0.5)

  f <- function(th) {
    M <- t(th[1] * X) # 3-by-2 transposed payload
    c(M[1, 1], M[2, 1], M[3, 1], M[1, 2], M[2, 2], M[3, 2])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # Each output is theta * t(X)[i, j]; the derivative wrt theta is t(X)[i, j],
  # flattened column-major.
  analytic <- matrix(as.vector(t(X)), ncol = 1)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("single-index subscripts a matrix-payload Ops result under exact autodiff", {
  # (theta[1] * X)[k] indexes column-major into the matrix payload rather than
  # raising a subscript-out-of-bounds error.
  X <- matrix(c(1, 2, 3, 4), nrow = 2)
  theta <- c(0.5)

  f <- function(th) {
    M <- th[1] * X
    c(M[1], M[2], M[3], M[4])
  }

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- matrix(as.vector(X), ncol = 1)
  expect_equal(exact, approx, tolerance = 1e-5)
  expect_equal(exact, analytic, tolerance = 1e-8)
})

test_that("length, dim, indexing, and t() work on a matrix-payload PrimalTangent", {
  pt <- primal_tangent(
    matrix(c(1, 2, 3, 4), nrow = 2),
    matrix(c(5, 6, 7, 8), nrow = 2)
  )
  expect_equal(length(pt), 4)
  expect_equal(dim(pt), c(2L, 2L))

  el <- pt[2]
  expect_s3_class(el, "PrimalTangent")
  expect_equal(el$primal, 2)
  expect_equal(el$tangent, 6)

  el_ij <- pt[1, 2]
  expect_s3_class(el_ij, "PrimalTangent")
  expect_equal(el_ij$primal, 3)
  expect_equal(el_ij$tangent, 7)

  tr <- t(pt)
  expect_s3_class(tr, "PrimalTangentArray")
  expect_equal(tr$primal, t(matrix(c(1, 2, 3, 4), nrow = 2)))
  expect_equal(tr$tangent, t(matrix(c(5, 6, 7, 8), nrow = 2)))
})

test_that("a vector-primal PrimalTangent recycles a broadcast tangent when indexed", {
  # A scalar broadcast tangent must expand to the primal length before selection.
  pt <- primal_tangent(c(10, 20, 30), 0.5)
  el <- pt[3]
  expect_s3_class(el, "PrimalTangent")
  expect_equal(el$primal, 30)
  expect_equal(el$tangent, 0.5)
})

# ---- Summary group methods on whole vectors and arrays ------------
#
# A Summary reduction (sum, prod, max, min) applied to the whole parameter vector
# (a PrimalTangentVector) or to a tangent-carrying array (a PrimalTangentArray,
# for example X %*% theta) must differentiate under exact mode. Only the scalar
# Summary.PrimalTangent method existed, so `sum(theta)` fell through to base sum
# on a list and errored with "invalid 'type' (list) of argument".
#
# Cross-checks follow the established autodiff convention: capprox at tolerance
# 1e-5 and the hand-derived analytic Jacobian at 1e-8.

test_that("autodiff differentiates sum over the whole parameter vector", {
  theta <- c(2, 5)
  f <- function(x) sum(x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # d sum(theta) / d theta_j = 1.
  expect_equal(exact, matrix(c(1, 1), nrow = 1), tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("autodiff differentiates sum over a scaled parameter vector", {
  theta <- c(2, 5)
  f <- function(x) sum(x * c(3, 4))

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # d sum(theta * c(3, 4)) / d theta_j = c(3, 4)[j].
  expect_equal(exact, matrix(c(3, 4), nrow = 1), tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("autodiff differentiates sum over a tangent array X %*% theta", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3) # 3-by-2 constant (data)
  theta <- c(0.5, -1)
  f <- function(x) sum(X %*% x)

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  # d sum(X theta) / d theta_j = colSums(X)[j].
  expect_equal(exact, matrix(colSums(X), nrow = 1), tolerance = 1e-8)
  expect_equal(exact, approx, tolerance = 1e-5)
})

test_that("Summary methods reduce a PrimalTangentVector", {
  ptv <- primal_tangent_vector(list(
    primal_tangent(2, 1),
    primal_tangent(5, 1),
    primal_tangent(3, 1)
  ))

  s <- sum(ptv)
  expect_s3_class(s, "PrimalTangent")
  expect_equal(s$primal, 10)
  expect_equal(s$tangent, 3)

  p <- prod(ptv)
  expect_equal(p$primal, 30)
  expect_equal(p$tangent, 5 * 3 + 2 * 3 + 2 * 5) # prod(x[-i]) summed over i

  mx <- max(ptv)
  expect_equal(mx$primal, 5)
  expect_equal(mx$tangent, 1)

  mn <- min(ptv)
  expect_equal(mn$primal, 2)
  expect_equal(mn$tangent, 1)
})

test_that("Summary methods reduce a PrimalTangentArray", {
  pta <- primal_tangent_array(c(2, 5, 3), c(1, 0, 0))

  s <- sum(pta)
  expect_s3_class(s, "PrimalTangent")
  expect_equal(s$primal, 10)
  expect_equal(s$tangent, 1)

  mx <- max(pta)
  expect_equal(mx$primal, 5)
  expect_equal(mx$tangent, 0)

  mn <- min(pta)
  expect_equal(mn$primal, 2)
  expect_equal(mn$tangent, 1)
})

test_that("Summary methods honor na.rm for prod, max, and min", {
  # A tangent container holding an NA primal element. With na.rm = TRUE the NA
  # element is dropped from both the primal and the tangent, matching base
  # Summary semantics on the primal while keeping the tangent aligned.
  pta <- primal_tangent_array(c(2, NA, 3), c(1, 9, 1))

  p <- prod(pta, na.rm = TRUE)
  expect_equal(p$primal, 6) # prod(c(2, 3))
  expect_equal(p$tangent, 3 * 1 + 2 * 1) # sum of prod(x[-i]) * t_i over kept i

  mx <- max(pta, na.rm = TRUE)
  expect_equal(mx$primal, 3)
  expect_equal(mx$tangent, 1)

  mn <- min(pta, na.rm = TRUE)
  expect_equal(mn$primal, 2)
  expect_equal(mn$tangent, 1)

  # sum stays consistent: the NA element is dropped from both slots.
  s <- sum(pta, na.rm = TRUE)
  expect_equal(s$primal, 5)
  expect_equal(s$tangent, 2)
})

# ---- tangent-stripping guard in extract_tangent_column ----

test_that("delta_method aborts when the transform strips tangents", {
  # base::rbind has no PrimalTangent method, so it builds a bare list matrix and
  # discards the tangents. This is the resolution any installed-package user gets
  # for a bare rbind() call from the global environment. The exact pass must
  # abort rather than return a silent all-zero Jacobian.
  strip <- function(theta) base::rbind(theta[1], theta[2])
  expect_error(
    delta_method(
      c(1, 2),
      transform = strip,
      covariance = diag(2),
      deriv_method = "exact"
    )
  )
})

test_that("delta_method with a tangent-preserving transform returns the derivative", {
  keep <- function(theta) c(theta[1], theta[2])
  result <- delta_method(
    c(1, 2),
    transform = keep,
    covariance = diag(2),
    deriv_method = "exact"
  )
  expect_equal(result, diag(2))
})

# ---- plain-list branch in pt_arrays ----

test_that("masked rbind of a c() pair list keeps tangents under exact mode", {
  # Inside the package namespace, c() on scalar pairs yields a plain list of
  # pairs; the masked rbind must route that through pt_flatten so the tangents
  # survive rather than collapsing to zero.
  f <- function(theta) rbind(c(2 * theta[1], 3 * theta[2]))
  result <- auto_differentiation(c(1, 1), f)
  expect_equal(result, diag(c(2, 3)))
})

# ---- out-of-bounds rejection in [.PrimalTangentVector ----

test_that("[.PrimalTangentVector rejects an out-of-bounds index", {
  v <- primal_tangent_vector(list(
    primal_tangent(1, 1),
    primal_tangent(2, 0)
  ))
  expect_error(v[3], "subscript out of bounds")
})

test_that("[.PrimalTangentVector still allows a negative index", {
  v <- primal_tangent_vector(list(
    primal_tangent(1, 1),
    primal_tangent(2, 0),
    primal_tangent(3, 0)
  ))
  dropped <- v[-3]
  expect_s3_class(dropped, "PrimalTangentVector")
  expect_equal(length(dropped), 2)
})

test_that("auto_differentiation aborts on an out-of-bounds parameter index", {
  f <- function(theta) theta[3]
  expect_error(auto_differentiation(c(1, 2), f), "subscript out of bounds")
})

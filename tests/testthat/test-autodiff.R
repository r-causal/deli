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
# Both argument forms are exercised. The multiple-argument form
# `prod(x[1], x[2], ...)` dispatches to `Summary.PrimalTangent`, and
# `prod(c(x[1], x[2], ...))` reduces a `PrimalTangentArray` through
# `Summary.PrimalTangentArray`, because `c()` on tangent-carrying values returns
# a tangent array rather than an unclassed list.
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

test_that("autodiff differentiates prod over a c() of parameters", {
  # `c()` on tangent-carrying values returns a PrimalTangentArray, so the
  # single-argument reduction dispatches to Summary.PrimalTangentArray and the
  # product rule applies across the concatenated elements.
  theta <- c(2, -3, 4)
  f <- function(x) prod(c(x[1], x[2], x[3]))

  exact <- auto_differentiation(theta, f)
  approx <- approx_differentiation(f, theta, method = "capprox")

  analytic <- rbind(c(-3 * 4, 2 * 4, 2 * -3)) # c(-12, 8, -6)
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

# ---- unsupported members on scalar pairs name the class ----------------------
#
# The scalar surfaces (`Ops.PrimalTangent`, `Math.PrimalTangent`, and the
# Summary reducer) route an unwired member to a cli abort that reports the
# PrimalTangent class rather than the array wording. These pin that rendering.

test_that("an unsupported binary operator on a scalar pair names the operator", {
  x <- primal_tangent(2, 1)
  expect_error(
    x %% x,
    'Operator "%%" is not supported for a <PrimalTangent>.',
    fixed = TRUE
  )
})

test_that("an unsupported unary operator on a scalar pair names the operator", {
  x <- primal_tangent(2, 1)
  expect_error(
    !x,
    'Unary "!" is not supported for a <PrimalTangent>.',
    fixed = TRUE
  )
})

test_that("an unsupported Math member on a scalar pair names the function", {
  x <- primal_tangent(2, 1)
  expect_error(
    gamma(x),
    'Math function "gamma" is not supported for a <PrimalTangent>.',
    fixed = TRUE
  )
})

test_that("an unsupported Summary member on a scalar pair names the function", {
  x <- primal_tangent(2, 1)
  expect_error(
    range(x),
    'Summary function "range" is not supported for a <PrimalTangent>.',
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

test_that("max and min on an all-NA container under na.rm return the base sentinel", {
  # With every primal element NA and na.rm = TRUE, na.rm drops the whole primal
  # to length zero. Base R reduces max(numeric(0)) to -Inf and min(numeric(0)) to
  # +Inf, each with a warning; the tangent surface must match that value and
  # wording. The sentinel is a constant, so its derivative is zero.
  pta <- primal_tangent_array(c(NA, NA), c(1, 1))

  expect_warning(
    mx <- max(pta, na.rm = TRUE),
    "no non-missing arguments to max; returning -Inf",
    fixed = TRUE
  )
  expect_equal(mx$primal, -Inf)
  expect_equal(mx$tangent, 0)

  expect_warning(
    mn <- min(pta, na.rm = TRUE),
    "no non-missing arguments to min; returning Inf",
    fixed = TRUE
  )
  expect_equal(mn$primal, Inf)
  expect_equal(mn$tangent, 0)
})

test_that("sum and prod on an all-NA container under na.rm return the empty reductions", {
  # An empty primal has well-defined sum (0) and product (1) via the base
  # reductions, each with a zero tangent and no warning.
  pta <- primal_tangent_array(c(NA, NA), c(1, 1))

  s <- sum(pta, na.rm = TRUE)
  expect_equal(s$primal, 0)
  expect_equal(s$tangent, 0)

  p <- prod(pta, na.rm = TRUE)
  expect_equal(p$primal, 1)
  expect_equal(p$tangent, 0)
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

test_that("the tangent-stripping abort recommends operations that carry tangents", {
  # The hint may name only operations that carry tangents from any environment,
  # which means the registered S3 methods: arithmetic, %*%, t(), and indexing.
  # matrix(), rbind(), cbind(), and ifelse() are masked inside the namespace and
  # so are unavailable to a user working in the global environment, which is why
  # the remedy is the t(X * resid) idiom and indicator arithmetic in place of
  # ifelse(). cli wraps the bullets, so the message is compared with runs of
  # whitespace collapsed to a single space.
  strip <- function(theta) base::rbind(theta[1], theta[2])
  msg <- tryCatch(
    auto_differentiation(c(1, 2), strip),
    error = conditionMessage
  )
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "carries no tangents")
  expect_match(flat, "t(X * resid)", fixed = TRUE)
  expect_match(flat, "ind * yes + (1 - ind) * no", fixed = TRUE)
  expect_match(flat, "capprox", fixed = TRUE)
})

test_that("the tangent-stripping abort does not recommend c()", {
  # A previous version of the hint advised combining tangent-carrying values
  # with c() as the remedy for a stripped result. c() does carry tangents, but
  # it cannot rebuild a p-by-n return, so recommending it here steered a user
  # from one failure into another.
  strip <- function(theta) base::rbind(theta[1], theta[2])
  msg <- tryCatch(
    auto_differentiation(c(1, 2), strip),
    error = conditionMessage
  )
  expect_false(grepl("`c()`", msg, fixed = TRUE))
})

test_that("a base ifelse in a transform aborts with the indicator-arithmetic advice", {
  # base::ifelse assigns each branch into a copy of `test`, which turns a
  # tangent-carrying branch into a bare list and discards the derivatives. base::
  # is named explicitly rather than forcing the transform's environment to the
  # global environment: pkgload attaches the masked internals to the search path
  # during devtools::test(), so an environment-only version of this test would
  # reach deli's tangent-aware ifelse, succeed, and assert nothing.
  branch <- function(theta) {
    base::ifelse(c(TRUE, FALSE, TRUE), theta[1] * 2, theta[2] * 3)
  }
  msg <- tryCatch(
    auto_differentiation(c(1, 2), branch),
    error = conditionMessage
  )
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "carries no tangents")
  expect_match(flat, "ind * yes + (1 - ind) * no", fixed = TRUE)
})

# ---- tangent-aware row-name assignment --------------------------------------

test_that("row-name assignment records the labels on a shaped tangent container", {
  # base::`rownames<-` is named explicitly to pin the route a user's estimating
  # function takes: it is base R's closure that reaches deli's `dimnames<-`
  # methods, and the test environment is a clone of the package namespace, so a
  # mask on the setter would otherwise stand in for the dispatch under test.
  #
  # Both slots carry the labels, not the primal alone. The subset methods index
  # the two slots with the same subscript, so a character subscript the primal
  # resolves and the tangent does not raises `no 'dimnames' attribute for array`
  # from the tangent. Labeling both slots is also what lets `dimnames()` read
  # the labels back, which is how a second assignment on the other axis keeps
  # the first.
  containers <- list(
    primal_tangent(
      base::matrix(c(1, 2, 3, 4), nrow = 2),
      base::matrix(0, 2, 2)
    ),
    primal_tangent_array(
      base::matrix(c(1, 2, 3, 4), nrow = 2),
      base::matrix(0, 2, 2)
    )
  )

  for (x in containers) {
    named <- base::`rownames<-`(x, c("a", "b"))
    expect_s3_class(named, class(x))
    expect_equal(dimnames(named), list(c("a", "b"), NULL))
    expect_equal(dimnames(named$primal), list(c("a", "b"), NULL))
    expect_equal(dimnames(named$tangent), list(c("a", "b"), NULL))
    # Labeling changes no value in either slot.
    expect_equal(as.vector(named$primal), as.vector(x$primal))
    expect_equal(as.vector(named$tangent), as.vector(x$tangent))
    # Assigning NULL to an unlabeled container returns within base R's own
    # frame, before either setter is reached, which is why the
    # `rownames(out) <- NULL` lines in R/ee-glm.R survive the exact pass.
    expect_identical(base::`rownames<-`(x, NULL), x)
    # Clearing labels that were recorded drops them from both slots. base R
    # leaves the all-NULL list its own setter builds rather than removing the
    # attribute, which is what a plain matrix reports too, so the container
    # answers no row labels exactly as the numeric pass does.
    cleared <- base::`rownames<-`(named, NULL)
    expect_equal(dimnames(cleared$primal), list(NULL, NULL))
    expect_equal(dimnames(cleared$tangent), list(NULL, NULL))
    expect_null(rownames(cleared))
  }
})

test_that("row-name assignment reshapes a broadcast tangent before labeling it", {
  # `theta[k] * X` carries a matrix primal against a scalar tangent, which has
  # no dimensions to label. That tangent broadcasts elementwise, so expanding it
  # to the primal's shape preserves every derivative it stands for and is what
  # lets the labels sit on both slots.
  x <- primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1)
  named <- base::`rownames<-`(x, c("a", "b"))
  expect_equal(
    named$tangent,
    base::matrix(1, 2, 2, dimnames = list(c("a", "b"), NULL))
  )
})

test_that("row-name assignment of the wrong length errors as base R does", {
  # Silencing the assignment silenced its validation with it: a label vector
  # longer than the axis is an error on the numeric pass, and the exact pass has
  # to stop in the same place rather than be the more permissive one.
  x <- primal_tangent_array(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    base::matrix(0, 2, 2)
  )
  expect_error(
    base::`rownames<-`(x, c("a", "b", "c")),
    "not equal to array extent"
  )
})

test_that("column names assigned after row names keep both axes", {
  # base R's `colnames<-` builds its replacement list from `dimnames(x)`, so the
  # labels already recorded have to read back or the second assignment erases
  # the first.
  x <- primal_tangent_array(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    base::matrix(0, 2, 2)
  )
  named <- base::`colnames<-`(
    base::`rownames<-`(x, c("a", "b")),
    c("u", "v")
  )
  expect_equal(dimnames(named), list(c("a", "b"), c("u", "v")))
  expect_equal(dimnames(named$primal), list(c("a", "b"), c("u", "v")))
  expect_equal(dimnames(named$tangent), list(c("a", "b"), c("u", "v")))
})

test_that("a labeled tangent container selects by name on both slots", {
  # The subset methods need no character-subscript rule of their own: base R
  # resolves the name against each slot's own dimnames, so recording the labels
  # on both slots is the whole of what name selection needs.
  x <- primal_tangent_array(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    base::matrix(c(10, 20, 30, 40), nrow = 2)
  )
  named <- base::`rownames<-`(x, c("a", "b"))
  row <- named["b", ]
  expect_equal(as.vector(row$primal), c(2, 4))
  expect_equal(as.vector(row$tangent), c(20, 40))

  pair <- base::`rownames<-`(
    primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1),
    c("a", "b")
  )
  expect_equal(as.vector(pair["a", ]$primal), c(1, 3))
  expect_equal(as.vector(pair["a", ]$tangent), c(1, 1))
})

test_that("a psi that labels its rows and selects by name agrees with the numeric pass", {
  # The two passes have to agree on this psi, and until the labels were recorded
  # they did not: the numeric pass resolved the names and the exact pass raised
  # `no 'dimnames' attribute for array` at the selection rather than at the
  # assignment. rbind() is the masked binder here, which is what an estimating
  # equation evaluated inside the package reaches, and it is the only route by
  # which a labeled container survives to be selected from at all.
  y <- c(1, 2, 3)
  psi <- function(theta) {
    out <- rbind(y - theta[1], (y - theta[1])^2 - theta[2])
    rownames(out) <- c("mean", "var")
    out[c("mean", "var"), ]
  }
  exact <- compute_bread(psi, c(2, 1), deriv_method = "exact")
  expect_equal(
    exact,
    compute_bread(psi, c(2, 1), deriv_method = "capprox"),
    tolerance = 1e-6
  )
  # The labels stay out of the Jacobian. The bread is read off the tangent
  # slots, flattened, so a labeled psi return reports the unnamed bread the
  # numeric pass reports.
  expect_null(dimnames(exact))
})

test_that("row-name assignment on a shapeless tangent container errors as base R does", {
  # base::`rownames<-` consults dim() before it reaches `dimnames<-`, so a
  # container that answers no dimensions never reaches deli's methods and
  # raises exactly the error a plain numeric vector raises. Both passes stop in
  # the same place rather than the exact pass being the more permissive one, so
  # nothing succeeds under exact differentiation that would fail without it.
  msg <- "attempt to set 'rownames' on an object with no dimensions"
  shapeless <- list(
    c(1, 2),
    primal_tangent(1, 2),
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  )

  for (x in shapeless) {
    expect_error(base::`rownames<-`(x, c("a", "b")), msg)
  }
})

# ---- plain-list branch in pt_arrays ----

test_that("masked rbind of a c() tangent array keeps tangents under exact mode", {
  # c() on scalar pairs yields a PrimalTangentArray, which the masked rbind
  # routes through pt_arrays' array branch, so the tangents survive the reshape.
  f <- function(theta) rbind(c(2 * theta[1], 3 * theta[2]))
  result <- auto_differentiation(c(1, 1), f)
  expect_equal(result, diag(c(2, 3)))
})

test_that("masked rbind of an lapply-produced pair list keeps tangents", {
  # pt_arrays keeps its plain-list branch for a list of scalar pairs built by
  # lapply() or sapply() inside a differentiated function; the masked rbind must
  # route that through pt_flatten so the tangents survive rather than collapsing
  # to zero.
  f <- function(theta) {
    rbind(lapply(seq_len(2), function(i) (i + 1) * theta[i]))
  }
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

# ---- c() on tangent-carrying values ----
#
# `c()` is a base internal generic, so an S3 method registered for it resolves
# from any environment, including a user's global environment. That makes it the
# one portable way to assemble a tangent-carrying vector: the masked reshaping
# helpers (matrix, rbind, cbind) are in scope only inside the package namespace.
# Every `c()` call that leads with a tangent-carrying operand returns a
# PrimalTangentArray whose slots are the concatenated primals and tangents, so
# the result feeds straight into the Ops, Math, and Summary methods.

test_that("c() on two scalar pairs concatenates both slots into a tangent array", {
  z <- c(primal_tangent(2, 1), primal_tangent(3, 0))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(2, 3))
  expect_equal(z$tangent, c(1, 0))
})

test_that("c() on two tangent arrays concatenates both slots", {
  z <- c(
    primal_tangent_array(c(1, 2), c(1, 0)),
    primal_tangent_array(c(3, 4), c(0, 1))
  )
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(1, 2, 3, 4))
  expect_equal(z$tangent, c(1, 0, 0, 1))
})

test_that("c() flattens a matrix-slotted tangent array column-major", {
  arr <- primal_tangent_array(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    base::matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2)
  )
  z <- c(arr)
  expect_s3_class(z, "PrimalTangentArray")
  expect_null(dim(z))
  expect_equal(z$primal, c(1, 2, 3, 4))
  expect_equal(z$tangent, c(0.1, 0.2, 0.3, 0.4))
})

test_that("c() on a PrimalTangentVector and a scalar pair keeps element order", {
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  z <- c(vec, primal_tangent(3, 0.5))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(1, 2, 3))
  expect_equal(z$tangent, c(1, 0, 0.5))
})

test_that("c() gives a plain numeric operand a zero tangent", {
  z <- c(primal_tangent(2, 1), 7, c(8, 9))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(2, 7, 8, 9))
  expect_equal(z$tangent, c(1, 0, 0, 0))
})

test_that("c() recycles a broadcast tangent to the primal length", {
  # A scalar pair built from a vector primal and a constant (the `theta[k] * X`
  # scaling idiom) carries a length-1 tangent that broadcasts under elementwise
  # arithmetic. Concatenation reshapes the two slots independently, so the
  # tangent is first expanded element for element.
  z <- c(primal_tangent(c(4, 5, 6), 2), primal_tangent(7, 1))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(4, 5, 6, 7))
  expect_equal(z$tangent, c(2, 2, 2, 1))
})

test_that("arithmetic on a c() result carries tangents", {
  z <- c(primal_tangent(2, 1), primal_tangent(3, 0)) * 2
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(4, 6))
  expect_equal(z$tangent, c(2, 0))
})

test_that("exp() on a c() result differentiates elementwise", {
  z <- exp(c(primal_tangent(0, 1), primal_tangent(1, 0)))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, exp(c(0, 1)))
  expect_equal(z$tangent, c(1, 0))
})

test_that("a c() result reaches auto_differentiation as a tangent array", {
  f <- function(x) c(x[1] * 2, x[2] * 3)
  expect_equal(auto_differentiation(c(1, 1), f), diag(c(2, 3)))
})

test_that("sum() over a c() result differentiates exactly", {
  f <- function(x) sum(c(x[1] * 2, x[2] * 3))
  expect_equal(auto_differentiation(c(1, 1), f), rbind(c(2, 3)))
})

test_that("exp() over a c() result differentiates exactly under a reduction", {
  theta <- c(0.3, -0.7)
  f <- function(x) sum(exp(c(x[1], x[2])))
  expect_equal(
    auto_differentiation(theta, f),
    rbind(exp(theta)),
    tolerance = 1e-12
  )
})

# ---- c()'s own formals are not data ----
#
# `c()` documents two formals, `recursive` and `use.names`. A method declared as
# `function(...)` swallows both into the dots and concatenates them as ordinary
# operands, which lengthens the result and, under differentiation, adds a
# spurious row to the Jacobian. The tangent methods declare both, as every base
# `c` method does.

test_that("c() matches use.names to its own formal rather than concatenating it", {
  z <- c(primal_tangent(2, 1), primal_tangent(3, 0), use.names = FALSE)
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(2, 3))
  expect_equal(z$tangent, c(1, 0))
})

test_that("c() matches recursive to its own formal rather than concatenating it", {
  z <- c(primal_tangent(2, 1), primal_tangent(3, 0), recursive = TRUE)
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(2, 3))
  expect_equal(z$tangent, c(1, 0))
})

test_that("a c() call passing use.names keeps the Jacobian shape", {
  f <- function(x) c(x[1] * 2, x[2] * 3, use.names = FALSE)
  expect_equal(auto_differentiation(c(1, 1), f), diag(c(2, 3)))
})

test_that("c() names the primal slot from its argument names and leaves the tangent bare", {
  # The tangent slot is a derivative array read positionally by
  # extract_tangent_column(), so it stays unnamed whatever the primal carries.
  z <- c(a = primal_tangent(2, 1), b = primal_tangent(3, 0))
  expect_equal(names(z$primal), c("a", "b"))
  expect_null(names(z$tangent))
})

test_that("c() drops the primal names when use.names is FALSE", {
  z <- c(a = primal_tangent(2, 1), b = primal_tangent(3, 0), use.names = FALSE)
  expect_null(names(z$primal))
  expect_equal(z$primal, c(2, 3))
})

# ---- c() dispatch limitations ----
#
# `c()` dispatches on its first argument only. Registering the three tangent
# classes therefore covers every call that leads with a tangent-carrying value,
# but a call that leads with a plain numeric reaches the internal default, which
# unpacks each pair into its two slots and returns a bare list. These tests pin
# that boundary so it is not later "fixed" by masking `c()` inside the namespace:
# a mask would only apply to code evaluated there and would leave a user's
# global-environment transform behaving differently.

test_that("c() leading with a plain numeric does not dispatch", {
  z <- c(1, primal_tangent(2, 1))
  expect_false(is_pt_array(z))
  expect_type(z, "list")
  # The internal default unpacked the pair into its primal and tangent slots
  expect_length(z, 3)
})

test_that("a numeric-first c() aborts rather than returning a zero Jacobian", {
  f <- function(theta) c(0, theta[1])
  expect_error(auto_differentiation(1, f), "carries no tangents")
})

test_that("c() dispatches when a tangent array or a tangent vector leads", {
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  vec <- primal_tangent_vector(list(primal_tangent(3, 0), primal_tangent(4, 1)))
  expect_s3_class(c(arr, 5), "PrimalTangentArray")
  expect_s3_class(c(vec, 5), "PrimalTangentArray")
})

# ---- coercing a tangent-carrying container to a plain double ----

test_that("as.numeric() on a tangent array aborts and names c() as the flatten", {
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  msg <- tryCatch(as.numeric(arr), error = conditionMessage)
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "tangent")
  expect_match(flat, "`c()`", fixed = TRUE)
})

test_that("as.numeric() on a PrimalTangentVector aborts with the same guidance", {
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  msg <- tryCatch(as.numeric(vec), error = conditionMessage)
  expect_match(gsub("[[:space:]]+", " ", msg), "tangent")
})

test_that("as.numeric() on a vector-payload scalar pair aborts", {
  # `theta[k] * X` scales a data vector by a single parameter, so both slots of
  # the scalar pair hold vectors. That value is as much a tangent-carrying
  # vector as a PrimalTangentArray is, and coercing it must abort rather than
  # hand back the primal with the derivative silently dropped.
  pair <- primal_tangent(c(2, 4, 6), 2)
  msg <- tryCatch(as.numeric(pair), error = conditionMessage)
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "tangent")
  expect_match(flat, "`c()`", fixed = TRUE)
})

test_that("as.numeric() on a matrix-payload scalar pair aborts", {
  pair <- primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1)
  expect_error(as.numeric(pair), "tangent-carrying")
})

test_that("coercing a scaled data vector aborts instead of zeroing the Jacobian", {
  f <- function(theta) sum(as.numeric(c(1, 2, 3) * theta[1]))
  expect_error(auto_differentiation(2, f), "tangent-carrying")
})

test_that("c() ports the as.numeric linear-predictor idiom to exact mode", {
  # The pkgdown articles flatten the n-by-1 matrix product with
  # inverse_logit(as.numeric(W %*% beta)). as.numeric() must return a double, so
  # it cannot carry a derivative and that idiom does not port to
  # deriv_method = "exact"; c() does, because it is a registered S3 method.
  #
  # What makes this a fair stand-in for a user's transform is the content of the
  # psi, not where it is evaluated: t(), %*%, c(), the arithmetic operators, and
  # inverse_logit() are none of them masked by deli, so each resolves to the
  # same function here as it would in a user's global environment. Scoping the
  # psi to globalenv() would establish nothing, because devtools::test() loads
  # with export_all = TRUE and so puts the masked helpers on the search path,
  # which globalenv()'s parent chain reaches.
  set.seed(4)
  n <- 60
  X <- cbind(1, rnorm(n), rbinom(n, 1, 0.4))
  beta <- c(-0.3, 0.8, 0.5)
  y <- rbinom(n, 1, inverse_logit(X %*% beta))

  psi <- function(theta) t(X * (y - inverse_logit(c(X %*% theta))))

  exact <- compute_bread(psi, beta, deriv_method = "exact")
  approx <- compute_bread(psi, beta, deriv_method = "capprox")
  expect_equal(exact, approx, tolerance = 1e-5)
})

# ---- mean() as a linear reduction -------------------------------------------

test_that("mean() on a PrimalTangentVector averages both slots", {
  vec <- primal_tangent_vector(list(
    primal_tangent(2, 1),
    primal_tangent(4, 0),
    primal_tangent(6, 0)
  ))
  result <- mean(vec)
  expect_s3_class(result, "PrimalTangent")
  expect_equal(result$primal, 4)
  expect_equal(result$tangent, 1 / 3)
})

test_that("mean() on a tangent array averages both slots", {
  arr <- primal_tangent_array(c(1, 2, 3, 4), c(1, 1, 0, 0))
  result <- mean(arr)
  expect_equal(result$primal, 2.5)
  expect_equal(result$tangent, 0.5)
})

test_that("mean() on a vector-payload scalar pair recycles the broadcast tangent", {
  # `theta[k] * X` leaves a scalar pair whose primal is a vector and whose
  # tangent is a single broadcast value. The mean has to average the tangent
  # over the primal's length rather than over the length-1 tangent slot, which
  # would divide a single derivative by the number of observations.
  pair <- primal_tangent(c(2, 4, 6), 2)
  result <- mean(pair)
  expect_equal(result$primal, 4)
  expect_equal(result$tangent, 2)
})

test_that("mean() raises no condition and agrees with sum(x) / length(x)", {
  # A tangent-carrying argument used to reach mean.default(), which warned
  # "argument is not numeric or logical: returning NA" and returned NA. The NA
  # then flowed into plain numeric arithmetic, so the Jacobian column read zero
  # and the standard error came back as zero.
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(3, 0)))
  expect_no_warning(mean(vec))
  by_hand <- sum(vec) / length(vec)
  expect_equal(mean(vec)$primal, by_hand$primal)
  expect_equal(mean(vec)$tangent, by_hand$tangent)
})

test_that("mean() with na.rm drops by primal NA and masks the tangent", {
  arr <- primal_tangent_array(c(1, NA, 3), c(1, 5, 0))
  result <- mean(arr, na.rm = TRUE)
  expect_equal(result$primal, 2)
  expect_equal(result$tangent, 0.5)
  expect_true(is.na(mean(arr)$primal))
})

test_that("mean() with a nonzero trim aborts", {
  arr <- primal_tangent_array(c(1, 2, 3, 4), c(1, 0, 0, 0))
  expect_error(
    mean(arr, trim = 0.25),
    class = "deli_exact_unsupported_function"
  )
})

test_that("mean() in a delta-method transform matches the finite difference", {
  covariance <- base::matrix(c(0.04, 0.01, 0.01, 0.09), 2, 2)
  transform <- function(theta) mean(theta)
  expect_equal(
    delta_method(
      c(1.5, 2.5),
      transform = transform,
      covariance = covariance,
      deriv_method = "exact"
    ),
    delta_method(
      c(1.5, 2.5),
      transform = transform,
      covariance = covariance,
      deriv_method = "capprox"
    ),
    tolerance = 1e-5
  )
})

test_that("mean() in a psi gives the same variance as the finite-difference bread", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) rbind(y - mean(theta))
  m <- m_estimate(stacked_equations = psi, init = 0, deriv_method = "exact")
  expect_equal(unname(coef(m)), mean(y))
  expect_gt(unname(sqrt(diag(vcov(m)))), 0)
  m_cap <- m_estimate(
    stacked_equations = psi,
    init = 0,
    deriv_method = "capprox"
  )
  expect_equal(vcov(m), vcov(m_cap), tolerance = 1e-6)
})

test_that("a base::rbind psi around mean() aborts instead of reporting a zero variance", {
  # deli's tangent-aware rbind() is masked inside the package namespace only,
  # and testthat::test_env() clones that namespace, so a bare rbind() in a test
  # resolves to a function a user in the global environment never reaches.
  # base:: is named explicitly to pin the resolution a user actually gets.
  #
  # Giving mean() a tangent rule turns this repro from a silent zero standard
  # error into a loud abort, which is an improvement but a different claim: a
  # base::rbind() of one tangent-carrying value builds a 1-by-2 list matrix, one
  # cell per slot, so the derivative is gone before the bread is summed.
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) base::rbind(y - mean(theta))
  expect_error(
    m_estimate(stacked_equations = psi, init = 0, deriv_method = "exact"),
    class = "deli_exact_tangent_lost"
  )
})

# ---- median() and quantile() decline order-statistic selection ---------------

test_that("median() aborts on every tangent surface", {
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(3, 0)))
  arr <- primal_tangent_array(c(1, 3), c(1, 0))
  pair <- primal_tangent(c(1, 3), 1)
  expect_error(median(vec), class = "deli_exact_unsupported_function")
  expect_error(median(arr), class = "deli_exact_unsupported_function")
  expect_error(median(pair), class = "deli_exact_unsupported_function")
})

test_that("quantile() aborts on every tangent surface", {
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(3, 0)))
  arr <- primal_tangent_array(c(1, 3), c(1, 0))
  pair <- primal_tangent(c(1, 3), 1)
  expect_error(quantile(vec, 0.5), class = "deli_exact_unsupported_function")
  expect_error(quantile(arr, 0.5), class = "deli_exact_unsupported_function")
  expect_error(quantile(pair, 0.5), class = "deli_exact_unsupported_function")
})

test_that("a median transform aborts rather than returning a 0-by-0 variance", {
  # median() on a tangent-carrying value used to return an empty
  # PrimalTangentVector, and delta_method handed back a 0-by-0 matrix with no
  # error and no NA to signal that anything had gone wrong.
  expect_error(
    delta_method(
      c(1.5, 2.5),
      transform = function(theta) median(theta),
      covariance = diag(2),
      deriv_method = "exact"
    ),
    class = "deli_exact_unsupported_function"
  )
})

test_that("a quantile transform aborts rather than differentiating an order statistic", {
  expect_error(
    delta_method(
      c(1.5, 2.5),
      transform = function(theta) quantile(theta, 0.5),
      covariance = diag(2),
      deriv_method = "exact"
    ),
    class = "deli_exact_unsupported_function"
  )
})

test_that("the order-statistic abort names the function and offers capprox", {
  arr <- primal_tangent_array(c(1, 3), c(1, 0))
  flat <- gsub(
    "[[:space:]]+",
    " ",
    tryCatch(median(arr), error = conditionMessage)
  )
  expect_match(flat, "`median()`", fixed = TRUE)
  expect_match(flat, "capprox", fixed = TRUE)
})

test_that("median() and quantile() on plain numeric data are unchanged", {
  y <- c(3, 1, 2, 5, 4)
  expect_equal(median(y), 3)
  expect_equal(unname(quantile(y, 0.25)), 2)
})

# ---- compiled base R functions name themselves under exact mode -------------

test_that("a compiled base R function names itself and its deli replacement", {
  err <- expect_error(
    auto_differentiation(0.5, function(theta) plogis(theta[1])),
    class = "deli_exact_unsupported_function"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "`plogis()`", fixed = TRUE)
  expect_match(flat, "inverse_logit", fixed = TRUE)
})

test_that("a namespace-qualified call is named without its namespace", {
  err <- expect_error(
    auto_differentiation(0.5, function(theta) stats::pnorm(theta[1])),
    class = "deli_exact_unsupported_function"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "`pnorm()`", fixed = TRUE)
  expect_match(flat, "standard_normal_cdf", fixed = TRUE)
})

test_that("a compiled function with no deli replacement is still named", {
  # Every distribution function raises the same compiled-code error, so the
  # abort reads the offender off the failing call and looks it up in a
  # replacement table rather than matching a fixed list of five names. qnorm()
  # has no deli counterpart and still has to be named, with the generic remedy.
  err <- expect_error(
    auto_differentiation(0.5, function(theta) qnorm(theta[1])),
    class = "deli_exact_unsupported_function"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "`qnorm()`", fixed = TRUE)
  expect_match(flat, "capprox", fixed = TRUE)
})

test_that("the psigamma abort records the reversed argument order", {
  err <- expect_error(
    auto_differentiation(1.5, function(theta) psigamma(theta[1], deriv = 1)),
    class = "deli_exact_unsupported_function"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "deli_polygamma", fixed = TRUE)
  expect_match(flat, "deriv", fixed = TRUE)
})

test_that("an unrelated error inside a differentiated function propagates unchanged", {
  boom <- function(theta) stop("a deliberate failure")
  err <- tryCatch(auto_differentiation(1, boom), error = function(e) e)
  expect_s3_class(err, "simpleError")
  expect_equal(conditionMessage(err), "a deliberate failure")
})

test_that("the tangent-loss abort keeps its class through the compiled-error handler", {
  strip <- function(theta) base::rbind(theta[1], theta[2])
  expect_error(
    auto_differentiation(c(1, 2), strip),
    class = "deli_exact_tangent_lost"
  )
})

# ---- coercing a genuinely scalar pair ---------------------------------------

test_that("as.numeric() on a genuinely scalar pair aborts", {
  # The last coercion hole. A scalar pair's primal is already a plain double, so
  # as.numeric() could hand it back and did, dropping the tangent with no NA and
  # no warning for any later rule to catch.
  pair <- primal_tangent(1.5, 1)
  expect_error(as.numeric(pair), class = "deli_exact_tangent_lost")
  expect_error(as.double(pair), class = "deli_exact_tangent_lost")
})

test_that("coercing a scalar pair aborts instead of zeroing the Jacobian", {
  transform <- function(theta) 2 * as.numeric(theta[1])
  expect_error(
    delta_method(
      c(1.5, 2.5),
      transform = transform,
      covariance = diag(2),
      deriv_method = "exact"
    ),
    class = "deli_exact_tangent_lost"
  )
  # The finite-difference reference is unaffected and returns the variance the
  # exact pass used to report as zero.
  expect_equal(
    delta_method(
      c(1.5, 2.5),
      transform = transform,
      covariance = diag(2),
      deriv_method = "capprox"
    ),
    base::matrix(4),
    tolerance = 1e-5
  )
})

test_that("as.logical() on a scalar pair still returns the primal logical", {
  # Coercion to logical stays, because a logical is not a type a derivative can
  # flow through: every value that reaches it has already dropped its tangent by
  # design, since comparisons return primal logicals and sign() has a zero
  # tangent rule. A double is different, because it is the type the tangent does
  # flow through, so returning one silently truncates a live derivative.
  expect_true(as.logical(primal_tangent(1, 1)))
  expect_false(as.logical(primal_tangent(0, 1)))
})

test_that("as.logical() reads the payload of every tangent-carrying container", {
  # Only the scalar pair had a method, so the two container surfaces fell
  # through to base R, which reads the classed list rather than the payload:
  # as.logical(primal_tangent_array(1, 1)) returned c(TRUE, TRUE), one value per
  # slot, and any longer payload raised `'list' object cannot be coerced to type
  # 'logical'`. The rule the scalar pair stands on covers all three, because a
  # logical coercion is a step function whose derivative is zero almost
  # everywhere.
  expect_identical(as.logical(primal_tangent_array(1, 1)), TRUE)
  expect_identical(
    as.logical(primal_tangent_array(c(1, 0, 2), c(1, 0, 1))),
    c(TRUE, FALSE, TRUE)
  )
  expect_identical(
    as.logical(primal_tangent_array(
      base::matrix(c(1, 0, 2, 0), nrow = 2),
      base::matrix(0, 2, 2)
    )),
    c(TRUE, FALSE, TRUE, FALSE)
  )
  expect_identical(
    as.logical(primal_tangent_vector(list(
      primal_tangent(1, 1),
      primal_tangent(0, 0),
      primal_tangent(2, 0)
    ))),
    c(TRUE, FALSE, TRUE)
  )
})

test_that("the indicator idiom differentiates with a container-valued condition", {
  # `ind * yes + (1 - ind) * no` is the conditional deli names in place of
  # ifelse(), and the condition an estimating function builds is a tangent
  # container whenever it came from theta, so the coercion has to answer for the
  # payload. Reading the slots instead returned a length-2 mask for a length-3
  # gate, which recycles rather than errors.
  f <- function(theta) {
    gate <- c(c(1, 0, 1) * theta[1])
    ind <- as.logical(gate)
    sum(ind * (theta[1] * c(1, 2, 3)) + (1 - ind) * (theta[2] * c(1, 2, 3)))
  }
  expect_equal(
    auto_differentiation(c(2, 3), f),
    base::matrix(c(4, 2), nrow = 1)
  )
  expect_equal(
    approx_differentiation(f, c(2, 3)),
    base::matrix(c(4, 2), nrow = 1),
    tolerance = 1e-5
  )
})

# ---- NA-scoped tangent-loss abort in extract_tangent_column ------------------

test_that("a plain numeric result containing NA aborts", {
  f <- function(theta) c(1, NA_real_)
  expect_error(auto_differentiation(1, f), class = "deli_exact_tangent_lost")
})

test_that("a tangent-free plain numeric result with no NA still reports zeros", {
  # A genuinely constant output has derivative zero and has to keep reporting
  # it. ee_percentile depends on the same rule: its indicator score carries no
  # tangent by construction, because comparisons return primal logicals by
  # design, and its identically zero bread is pinned in
  # test-exact-mode-acceptance.R.
  f <- function(x) 5
  expect_equal(auto_differentiation(1, f)[1, 1], 0)
})

test_that("a list matrix from base::rbind aborts rather than failing inside rowSums", {
  # base::rbind() on one tangent-carrying value builds a 1-by-2 list matrix, one
  # cell per slot. It carries a dim, so it misses compute_bread's list-of-pairs
  # branch and used to reach rowSums(), which failed with the opaque
  # "'x' must be numeric or complex" before any diagnostic could run.
  y <- c(1, 2, 3)
  psi <- function(theta) base::rbind(y - theta[1])
  expect_error(
    compute_bread(psi, 0, deriv_method = "exact"),
    class = "deli_exact_tangent_lost"
  )
})

# ---- a length-1 payload that carries a dim ----------------------------------

test_that("a length-1 payload with a dim is not treated as a scalar pair", {
  # pt_is_scalar() rejects a payload that carries dimensions even when it holds
  # a single value, which is the clause a naive length-1 test would drop. With
  # the dim clause the pair takes the array paths, so `[` subsets the payload
  # and returns a dimensionless scalar pair; without it the pair would
  # self-select and hand back a primal that is still a 1-by-1 matrix.
  pair <- primal_tangent(base::matrix(5, 1, 1), 1)
  expect_false(pt_is_scalar(pair))
  expect_equal(dim(pair), c(1L, 1L))
  selected <- pair[1]
  expect_s3_class(selected, "PrimalTangent")
  expect_null(dim(selected$primal))
  expect_equal(selected$primal, 5)
  expect_equal(selected$tangent, 1)
  expect_error(as.numeric(pair), class = "deli_exact_tangent_lost")
})

# ---- the compiled-error rewrite is scoped to functions that never dispatch ---

test_that("the Math group listing separates dispatching members from stats functions", {
  # Every entry in the replacement table has to stay outside the group, or its
  # rewrite would be skipped.
  expect_true(all(c("log", "sqrt", "round") %in% pt_math_group_members))
  expect_false(any(
    c("psigamma", "plogis", "qlogis", "pnorm", "dnorm", "qnorm") %in%
      pt_math_group_members
  ))
  # Every member listed dispatches, which is what makes the compiled-code
  # message from one of them a statement about its data rather than about deli.
  # A primal of 1 is inside the domain of every member, so none of them warns
  # about producing a NaN.
  pair <- primal_tangent(1, 1)
  dispatches <- vapply(
    pt_math_group_members,
    function(fname) {
      result <- tryCatch(do.call(fname, list(pair)), error = conditionMessage)
      is_pt(result) || grepl("is not supported for a", result)
    },
    logical(1)
  )
  expect_true(all(dispatches))
})

test_that("a Math group member failing on its data keeps the base error", {
  # log() dispatches, so a tangent-carrying argument reaches deli's own method
  # and never gets as far as compiled code. The same message raised from a
  # character column is a data problem, and rewriting it would claim deli
  # cannot differentiate a function it differentiates.
  dose <- c("low", "high")
  psi <- function(theta) t(log(dose) - theta[1])
  err <- tryCatch(
    compute_bread(psi, 0, deriv_method = "exact"),
    error = function(e) e
  )
  expect_false(inherits(err, "deli_exact_unsupported_function"))
  expect_match(
    conditionMessage(err),
    "non-numeric argument to mathematical function"
  )
})

test_that("a compiled function outside the group generics is still rewritten", {
  err <- expect_error(
    auto_differentiation(0.5, function(theta) plogis(theta[1])),
    class = "deli_exact_unsupported_function"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "`plogis()`", fixed = TRUE)
})

test_that("a compiled call reached indirectly keeps the base error", {
  # The failing call is `get("plogis")(theta[1])`, whose function position is a
  # call rather than a name, so no offender can be read off it. Naming the
  # wrong function would be worse than handing back what base R said.
  err <- tryCatch(
    auto_differentiation(0.5, function(theta) get("plogis")(theta[1])),
    error = function(e) e
  )
  expect_false(inherits(err, "deli_exact_unsupported_function"))
  expect_match(
    conditionMessage(err),
    "[Nn]on-numeric argument to mathematical function"
  )
})

# ---- NaN is a constant, not evidence of a dropped tangent -------------------

test_that("a tangent-free NaN reports zero derivatives rather than aborting", {
  # anyNA() counts NaN as missing, which would report a lost derivative for a
  # value arithmetic on numbers produces. The NA rule is about a function that
  # returned NA for an argument it did not recognize.
  f <- function(theta) c(1, NaN)
  expect_no_error(auto_differentiation(1, f))
  expect_equal(auto_differentiation(1, f)[, 1], c(0, 0))
})

# ---- mean() over a matrix payload -------------------------------------------

test_that("mean() on a matrix-payload tangent array averages every cell", {
  arr <- primal_tangent_array(
    base::matrix(c(1, 2, 3, 4), 2, 2),
    base::matrix(c(1, 1, 0, 0), 2, 2)
  )
  result <- mean(arr)
  expect_s3_class(result, "PrimalTangent")
  expect_equal(result$primal, 2.5)
  expect_equal(result$tangent, 0.5)
  expect_null(dim(result$primal))
})

test_that("mean() on a tangent array matches base mean() on the primal", {
  # sum(x) / n and mean(x) disagree around the tenth significant digit on a
  # cancellation-heavy vector, because mean() makes a second pass to correct the
  # rounding error. The primal slot follows base R.
  values <- c(1e16, 1, -1e16, 2, 3)
  arr <- primal_tangent_array(values, rep(1, 5))
  expect_identical(mean(arr)$primal, mean(values))
})

# ---- intact tangents in an unsupported container ----------------------------

test_that("a list of tangent arrays reports an unsupported shape", {
  # An lapply() over several equations returns a list whose elements each still
  # carry their tangents. Nothing was lost, so reporting a lost tangent would
  # send the reader looking for a stripping function that is not there.
  X <- base::matrix(c(1, 1, 1, 2), 2, 2)
  psi <- function(theta) lapply(1:2, function(k) (X %*% theta) * k)
  err <- expect_error(
    compute_bread(psi, c(0, 0), deriv_method = "exact"),
    class = "deli_exact_unsupported_shape"
  )
  flat <- gsub("[[:space:]]+", " ", conditionMessage(err))
  expect_match(flat, "container shape", fixed = TRUE)
})

test_that("a list matrix with no tangents still reports a lost tangent", {
  # base::rbind() on one tangent-carrying value builds a list matrix of plain
  # numeric slots, so the derivative really is gone and the two conditions stay
  # distinguishable.
  y <- c(1, 2, 3)
  psi <- function(theta) base::rbind(y - theta[1])
  expect_error(
    compute_bread(psi, 0, deriv_method = "exact"),
    class = "deli_exact_tangent_lost"
  )
})

test_that("a list-shaped return under capprox is left alone", {
  # The list guard is scoped to the exact pass. A data.frame return is
  # list-shaped and carries a dim, so rowSums() reduces it as intended under
  # finite differences, and the guard must not intercept it.
  y <- c(1, 2, 3, 4)
  psi <- function(theta) as.data.frame(base::rbind(y - theta[1]))
  expect_equal(
    compute_bread(psi, 2, deriv_method = "capprox"),
    base::matrix(4),
    tolerance = 1e-6
  )
})

# ---- coercing a tangent-carrying container to a plain matrix -----------------

test_that("as.matrix() aborts on every tangent-carrying container", {
  # A matrix of doubles is a type the tangent flows through, so a coercion to one
  # can only truncate a live derivative. That is the ground as.double() stands
  # on, and the same abort covers this coercion.
  #
  # What base R does instead is decided by shapes that have nothing to do with
  # the payload: as.matrix.default() reads length() and dim() through deli's
  # methods, which answer for the payload, while is.matrix() reads the real
  # attribute and sees a classed list of two slots. The two disagree, so the
  # coercion raises `dims [product N] do not match the length of object [2]` from
  # a base frame no deli guard runs in.
  containers <- list(
    primal_tangent(2.7, 1),
    primal_tangent(c(2, 4, 6), 2),
    primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1),
    primal_tangent_array(
      base::matrix(c(1, 2, 3, 4), nrow = 2),
      base::matrix(0, 2, 2)
    ),
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  )

  for (x in containers) {
    expect_error(as.matrix(x), class = "deli_exact_tangent_lost")
  }
})

test_that("as.matrix() aborts on the payloads base R silently reshapes", {
  # A two-element payload is the one shape as.matrix.default() accepts, because
  # the container itself holds two slots: it sets `dim = c(2, 1)` on the pair and
  # hands back an object that still claims the tangent class while its cells are
  # the primal and tangent slots. A one-parameter PrimalTangentVector reshapes
  # the same way, holding one slot. Neither result is a matrix of anything, and
  # nothing downstream reports that.
  expect_error(
    as.matrix(primal_tangent_array(c(1, 2), c(1, 0))),
    class = "deli_exact_tangent_lost"
  )
  expect_error(
    as.matrix(primal_tangent(c(1, 2), c(1, 0))),
    class = "deli_exact_tangent_lost"
  )
  expect_error(
    as.matrix(primal_tangent_vector(list(primal_tangent(1, 1)))),
    class = "deli_exact_tangent_lost"
  )
})

test_that("the as.matrix() abort names what to do instead", {
  # The abort follows the shape of the one as.double() raises: it says what was
  # lost and what to reach for instead, so a reader who wrote the coercion has a
  # route out of it rather than only a refusal.
  msg <- tryCatch(as.matrix(primal_tangent(2.7, 1)), error = conditionMessage)
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "tangent")
  expect_match(flat, "needs no coercion")
})

test_that("as.matrix() in a psi aborts rather than zeroing the Jacobian", {
  # The two-element payload is the dangerous one: base R reshapes it, the
  # tangents are gone from what comes back, and the exact pass reports zero where
  # finite differences report 3.
  f <- function(theta) sum(as.matrix(c(1, 2) * theta[1]))
  expect_error(auto_differentiation(2, f), class = "deli_exact_tangent_lost")
  expect_equal(
    approx_differentiation(f, 2),
    base::matrix(3),
    tolerance = 1e-5
  )
})

test_that("as.matrix() in a psi reports the tangent loss, not the dims message", {
  y <- c(1, 2, 3)
  psi <- function(theta) t(as.matrix(y - theta[1]))
  err <- expect_error(
    compute_bread(psi, 1, deriv_method = "exact"),
    class = "deli_exact_tangent_lost"
  )
  expect_false(grepl("do not match the length", conditionMessage(err)))
})

test_that("as.matrix() on ordinary values is untouched", {
  expect_equal(as.matrix(c(1, 2, 3)), base::matrix(c(1, 2, 3), 3, 1))
  m <- base::matrix(c(1, 2, 3, 4), nrow = 2)
  expect_identical(as.matrix(m), m)
  expect_equal(
    as.matrix(data.frame(a = c(1, 2))),
    base::matrix(c(1, 2), 2, 1, dimnames = list(NULL, "a"))
  )
})

# ---- coercing a tangent-carrying container to an integer ---------------------

test_that("as.integer() on a scalar pair aborts instead of returning both slots", {
  # The corruption this closes: base R coerced the classed list itself, so
  # as.integer(primal_tangent(2.7, 1)) returned c(2L, 1L), the truncated primal
  # with the tangent appended to it as a second value. Nothing downstream can
  # tell that from data.
  expect_error(
    as.integer(primal_tangent(2.7, 1)),
    class = "deli_exact_tangent_lost"
  )
})

test_that("as.integer() aborts on the shapes base R refused to coerce", {
  # These already stopped, with `'list' object cannot be coerced to type
  # 'integer'` from a base frame, which points at neither the cause nor the
  # remedy. The abort replaces it so every shape reports the same rule.
  containers <- list(
    primal_tangent(c(2, 4, 6), 2),
    primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1),
    primal_tangent_array(c(1, 2), c(1, 0)),
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  )

  for (x in containers) {
    expect_error(as.integer(x), class = "deli_exact_tangent_lost")
  }
})

test_that("the as.integer() abort names what to do instead", {
  # Same guidance as the coercion to a double, for the same reason: an integer is
  # a plain number, so the tangent has nowhere to go.
  msg <- tryCatch(
    as.integer(primal_tangent(2.7, 1)),
    error = conditionMessage
  )
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "tangent")
  expect_match(flat, "needs no coercion")
})

test_that("as.integer() in a psi aborts rather than zeroing the Jacobian", {
  # The silent coercion took the derivative with it: the exact pass reported zero
  # where finite differences report the truncation-free value.
  f <- function(theta) sum(as.integer(theta[1]) * 2)
  expect_error(auto_differentiation(2.7, f), class = "deli_exact_tangent_lost")
})

test_that("as.vector() with a numeric mode aborts on a tangent-carrying value", {
  # as.vector() dispatches on the class of its first argument, so the rule
  # reaches the mode spellings of the coercion too. Each of these was a way
  # around a method that already aborts: as.vector(x, "integer") returned
  # c(2L, 1L) exactly as as.integer() did, and as.vector(x, "double") returned
  # c(2.7, 1) even though as.double() itself stops.
  pair <- primal_tangent(2.7, 1)
  for (target in c("integer", "double", "numeric")) {
    expect_error(
      as.vector(pair, target),
      class = "deli_exact_tangent_lost"
    )
  }
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_error(as.vector(arr, "integer"), class = "deli_exact_tangent_lost")
  expect_error(as.vector(arr, "double"), class = "deli_exact_tangent_lost")
})

test_that("as.vector() with no mode hands the container back with its tangents", {
  # The default mode is "any", which on a list returns the list, so
  # `as.vector(X %*% theta)` keeps its derivatives and differentiates exactly.
  # That has to survive: the abort is about the modes that force an atomic type.
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_identical(as.vector(arr), arr)
  f <- function(theta) sum(as.vector(c(1, 2, 3) * theta[1]))
  expect_equal(auto_differentiation(2, f)[1, 1], 6)
})

test_that("integer coercion of ordinary values is untouched", {
  expect_identical(as.integer(c(1.7, 2.2)), c(1L, 2L))
  expect_identical(as.vector(c(1.7, 2.2), "integer"), c(1L, 2L))
  expect_identical(as.vector(c(1.7, 2.2), "double"), c(1.7, 2.2))
})

# ---- coercing a tangent-carrying container to a character or a complex -------

test_that("as.character() aborts instead of rendering the two slots", {
  # The same defect family as the integer coercion, in the type whose corruption
  # is hardest to notice: as.character(primal_tangent(2.7, 1)) returned
  # c("2.7", "1"), the primal with the tangent appended to it as a second value,
  # and on a container it deparsed each slot into a string of its own, so
  # "c(1, 2)" stood where a value was expected.
  containers <- list(
    primal_tangent(2.7, 1),
    primal_tangent(c(2, 4, 6), 2),
    primal_tangent_array(c(1, 2), c(1, 0)),
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  )

  for (x in containers) {
    expect_error(as.character(x), class = "deli_exact_tangent_lost")
  }
})

test_that("as.complex() aborts on a tangent-carrying container", {
  # A complex has a real and an imaginary part, neither of which is a place to
  # keep a derivative, so it is a plain number for this purpose. The scalar pair
  # returned c(2.7+0i, 1+0i); the wider payloads already stopped inside base R.
  containers <- list(
    primal_tangent(2.7, 1),
    primal_tangent_array(c(1, 2), c(1, 0)),
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  )

  for (x in containers) {
    expect_error(as.complex(x), class = "deli_exact_tangent_lost")
  }
})

test_that("as.vector() aborts for the character and complex modes too", {
  # Each mode spelling is another way of asking for the coercion the method
  # above refuses, and each was a way around it: as.vector(x, "character")
  # returned c("2.7", "1") exactly as as.character() did.
  pair <- primal_tangent(2.7, 1)
  for (target in c("character", "complex")) {
    expect_error(as.vector(pair, target), class = "deli_exact_tangent_lost")
  }
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_error(as.vector(arr, "character"), class = "deli_exact_tangent_lost")
  expect_error(as.vector(arr, "complex"), class = "deli_exact_tangent_lost")
})

test_that("as.vector() with mode 'list' still hands the slots back", {
  # A list is the one non-numeric mode a container has a faithful representation
  # in, since it is what the container already is, so this mode keeps
  # delegating to base R and returns both slots with no coercion at all.
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_identical(
    as.vector(arr, "list"),
    list(primal = c(1, 2), tangent = c(1, 0))
  )
})

test_that("character and complex coercion of ordinary values is untouched", {
  expect_identical(as.character(c(1.5, 2)), c("1.5", "2"))
  expect_identical(as.complex(c(1.5, 2)), c(1.5 + 0i, 2 + 0i))
  expect_identical(as.vector(c(1.5, 2), "character"), c("1.5", "2"))
  expect_identical(as.vector(c(1.5, 2), "complex"), c(1.5 + 0i, 2 + 0i))
})

test_that("a character coercion inside a psi aborts rather than zeroing the bread", {
  # The silent version reported a zero derivative for a parameter the function
  # genuinely depends on, because the string it built carried no tangent for
  # anything downstream to read.
  f <- function(theta) sum(as.numeric(as.character(theta[1])) * 2)
  expect_error(auto_differentiation(2.7, f), class = "deli_exact_tangent_lost")
})

# ---- the replacement functions that reach the coercion methods ---------------

test_that("mode<- reaches the coercion abort", {
  # `mode<-` is an ordinary closure that looks up `as.<value>` and calls it, so
  # every method above governs the assignment spelling of the same coercion
  # without a setter of its own.
  numeric_pair <- primal_tangent(2.7, 1)
  expect_error(
    mode(numeric_pair) <- "numeric",
    class = "deli_exact_tangent_lost"
  )
  integer_pair <- primal_tangent(2.7, 1)
  expect_error(
    mode(integer_pair) <- "integer",
    class = "deli_exact_tangent_lost"
  )
  character_pair <- primal_tangent(2.7, 1)
  expect_error(
    mode(character_pair) <- "character",
    class = "deli_exact_tangent_lost"
  )
})

test_that("storage.mode<- corrupts a scalar pair, which is the one hole no method closes", {
  # `storage.mode<-` is a primitive that coerces the object it is handed without
  # consulting the S3 method table: registering a method for it changes nothing,
  # and it does not route through `mode<-`, which is what does reach the abort.
  # So the corruption is pinned rather than fixed. The result is a classed
  # integer vector holding the truncated primal beside the tangent, which is the
  # same two-slots-as-data shape `as.integer()` used to return.
  x <- primal_tangent(2.7, 1)
  storage.mode(x) <- "integer"
  expect_s3_class(x, "PrimalTangent")
  expect_identical(unclass(x), c(primal = 2L, tangent = 1L))

  # The wider payloads are refused by base R itself, for want of a coercion from
  # a list rather than by any rule of deli's.
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_error(storage.mode(arr) <- "integer")
})

# ---- two-index selection when the tangent carries no dimensions --------------

test_that("the two-index form selects from a length-1 matrix payload", {
  # pt_recycle_tangent() reshaped only a tangent shorter than its primal, so a
  # 1-by-1 payload against a scalar tangent left the tangent dimensionless and
  # `tangent[i, j, drop = FALSE]` raised `incorrect number of dimensions`. The
  # primal's dim is what decides the shape of both slots, whatever their lengths.
  pair <- primal_tangent(base::matrix(5, 1, 1), 1)
  selected <- pair[1, 1]
  expect_s3_class(selected, "PrimalTangent")
  expect_null(dim(selected$primal))
  expect_equal(selected$primal, 5)
  expect_equal(selected$tangent, 1)
})

test_that("the two-index form selects from an equal-length flat tangent", {
  # The same clause covers every equal-length case and not only the length-1 one:
  # a 2-by-2 payload against a length-4 dimensionless tangent raised the same
  # error for `[i, j]`, for `[i, ]`, and for `[, j]`. Both slots are pinned, so a
  # tangent that is reshaped to the wrong layout is caught rather than passing
  # for a fix.
  pair <- primal_tangent(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    c(10, 20, 30, 40)
  )

  cell <- pair[1, 2]
  expect_s3_class(cell, "PrimalTangent")
  expect_equal(cell$primal, 3)
  expect_equal(cell$tangent, 30)

  row <- pair[2, ]
  expect_equal(row$primal, base::matrix(c(2, 4), 1, 2))
  expect_equal(row$tangent, base::matrix(c(20, 40), 1, 2))

  column <- pair[, 2]
  expect_equal(column$primal, base::matrix(c(3, 4), 2, 1))
  expect_equal(column$tangent, base::matrix(c(30, 40), 2, 1))
})

test_that("single-index selection and a broadcast tangent are unchanged", {
  broadcast <- primal_tangent(base::matrix(c(1, 2, 3, 4), nrow = 2), 1)
  cell <- broadcast[1, 2]
  expect_equal(cell$primal, 3)
  expect_equal(cell$tangent, 1)

  flat <- primal_tangent(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    c(10, 20, 30, 40)
  )
  first <- flat[1]
  expect_equal(first$primal, 1)
  expect_equal(first$tangent, 10)
})

test_that("a vector payload still rejects a second subscript", {
  # A payload with no dimensions keeps base R's answer to a two-index selection,
  # so the exact pass stops where a plain numeric vector stops.
  pair <- primal_tangent(c(2, 4, 6), 2)
  expect_equal(pair[2]$primal, 4)
  expect_equal(pair[2]$tangent, 2)
  expect_error(pair[2, 1], "incorrect number of dimensions")
})

test_that("two-index selection of a 1-by-1 payload keeps working in a psi", {
  # `X * theta[k]` on a 1-by-1 design is the route a psi reaches this shape by,
  # and the arithmetic gives both slots the design's dim, so this selection
  # already worked. It pins the natural path against a fix that reshapes the
  # wrong slot.
  design <- base::matrix(2, 1, 1)
  f <- function(theta) sum((design * theta[1])[1, 1])
  expect_equal(auto_differentiation(3, f)[1, 1], 2)
})

# ---- container recognition across the three tangent classes ------------------
#
# `is_tangent_container()` is the predicate the masked reshaping helpers, the
# masked binders, and `c()` consult before deciding whether an operand carries
# derivatives, so an operand it does not recognize takes the numeric-constant
# route and its tangents are replaced with zeros. A list is one of the shapes an
# `lapply()` inside a differentiated function produces, and which of the three
# tangent classes its elements hold is not the author's choice: `theta[k] * x`
# yields a scalar pair, `c(theta[k] * x)` yields a tangent array, and
# `theta[i:j]` yields a parameter vector. All three carry derivatives, so a list
# of any of them is a container.

test_that("a list of tangent arrays is a tangent container", {
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  expect_true(is_tangent_container(list(arr, arr)))
})

test_that("a list of parameter vectors is a tangent container", {
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  expect_true(is_tangent_container(list(vec, vec)))
})

test_that("a list mixing the tangent classes is a tangent container", {
  pair <- primal_tangent(2, 1)
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  expect_true(is_tangent_container(list(arr, vec)))
  expect_true(is_tangent_container(list(vec, arr)))
  # A scalar pair anywhere in the list was already recognized, whichever
  # position it took.
  expect_true(is_tangent_container(list(pair, arr)))
  expect_true(is_tangent_container(list(arr, pair)))
})

test_that("each tangent class is a container on its own", {
  expect_true(is_tangent_container(primal_tangent(2, 1)))
  expect_true(is_tangent_container(primal_tangent_array(c(1, 2), c(1, 0))))
  expect_true(is_tangent_container(
    primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  ))
})

test_that("a list holding no tangent-carrying element is not a container", {
  # A data frame is list-shaped and reaches this predicate on the
  # finite-difference pass, where it has to keep its numeric route: rowSums()
  # reduces it as intended, which is what pins a list-shaped return under
  # capprox above.
  expect_false(is_tangent_container(list(1, 2)))
  expect_false(is_tangent_container(list(a = 1, b = list(2, 3))))
  expect_false(is_tangent_container(list()))
  expect_false(is_tangent_container(data.frame(a = c(1, 2))))
})

test_that("a plain value is not a tangent container", {
  expect_false(is_tangent_container(c(1, 2)))
  expect_false(is_tangent_container(base::matrix(c(1, 2, 3, 4), nrow = 2)))
  expect_false(is_tangent_container(NULL))
  expect_false(is_tangent_container("a"))
})

# ---- a list operand c() cannot flatten ---------------------------------------
#
# Both slots of a tangent-carrying value are plain numeric objects of identical
# shape, and `pt_concat()` keeps that contract by normalizing every operand to
# parallel primal and tangent vectors of the same length. A list operand breaks
# it silently: `as.vector()` hands a list back unchanged, so the fallback pairs
# it with one zero per top-level element while the `unlist()` that follows
# flattens its internals into the primal slot. The two slots then have different
# lengths, and the Jacobian read off the tangent slot has fewer rows than the
# value it claims to differentiate, which the delta method goes on to multiply
# into a covariance matrix of the wrong size.
#
# The one list shape that does flatten is a list of scalar pairs, which is what
# an `lapply()` over `theta[k]` produces, and it keeps working. Any other list
# still holds whatever derivatives it was given, so it is a container shape with
# no reduction rather than a lost tangent, and the abort names
# `do.call(c, ...)` as the concatenation that does work.

test_that("c() aborts on a list operand rather than mis-shaping its slots", {
  expect_error(
    c(primal_tangent(2, 1), list(a = 1, b = list(2, 3))),
    class = "deli_exact_unsupported_shape"
  )
})

test_that("c() aborts on a list of tangent-carrying values it cannot flatten", {
  arr <- primal_tangent_array(c(1, 2), c(1, 0))
  vec <- primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
  expect_error(
    c(primal_tangent(2, 1), list(arr, arr)),
    class = "deli_exact_unsupported_shape"
  )
  expect_error(
    c(primal_tangent(2, 1), list(vec, vec)),
    class = "deli_exact_unsupported_shape"
  )
  expect_error(
    c(primal_tangent(2, 1), list(primal_tangent(3, 0), arr)),
    class = "deli_exact_unsupported_shape"
  )
})

test_that("the list-operand abort names the concatenation that does work", {
  # Reporting a lost tangent here would misdiagnose it twice over: nothing was
  # lost, and the remedy that abort gives is to reach for c(), which is the call
  # that failed. cli wraps the bullets, so the message is compared with runs of
  # whitespace collapsed to a single space.
  msg <- tryCatch(
    c(primal_tangent(2, 1), list(a = 1, b = list(2, 3))),
    error = conditionMessage
  )
  flat <- gsub("[[:space:]]+", " ", msg)
  expect_match(flat, "do.call(c, ...)", fixed = TRUE)
})

test_that("c() on a list of scalar pairs still concatenates both slots", {
  z <- c(primal_tangent(2, 1), list(primal_tangent(3, 0), primal_tangent(4, 1)))
  expect_s3_class(z, "PrimalTangentArray")
  expect_equal(z$primal, c(2, 3, 4))
  expect_equal(z$tangent, c(1, 0, 1))
})

test_that("a transform concatenating a list of groups aborts under exact mode", {
  # The shape a user reaches this by: one tangent-carrying value per group from
  # an lapply(), concatenated onto a parameter with c(). The tangent slot then
  # holds one zero per group where the primal holds every group's values, so the
  # Jacobian comes back with a handful of rows rather than one per output and the
  # delta method returns a covariance matrix of that size.
  xs <- list(c(1, 2), c(3, 4), c(5, 6))
  concat_groups <- function(theta) {
    per_group <- lapply(xs, function(x) c(theta[2] * x))
    c(theta[1], per_group)
  }
  expect_error(
    auto_differentiation(c(0.5, 2), concat_groups),
    class = "deli_exact_unsupported_shape"
  )
  expect_error(
    delta_method(
      c(0.5, 2),
      transform = concat_groups,
      covariance = diag(2),
      deriv_method = "exact"
    ),
    class = "deli_exact_unsupported_shape"
  )
  # The finite-difference pass refuses the same transform, because c() on a
  # plain numeric and a list returns a list rather than a vector. Nothing may
  # succeed under exact differentiation that fails without it.
  expect_error(
    delta_method(c(0.5, 2), transform = concat_groups, covariance = diag(2)),
    "cannot be coerced"
  )
})

test_that("do.call(c, ...) over a list of groups differentiates exactly", {
  # The supported route, which the abort above names. Every operand reaches
  # pt_concat() in its own right rather than inside a list, so both slots grow
  # together and the Jacobian has one row per output.
  xs <- list(c(1, 2), c(3, 4), c(5, 6))
  concat_groups <- function(theta) {
    per_group <- lapply(xs, function(x) c(theta[2] * x))
    do.call(c, c(list(theta[1]), per_group))
  }
  jacobian <- auto_differentiation(c(0.5, 2), concat_groups)
  expect_equal(dim(jacobian), c(7L, 2L))
  expect_equal(
    jacobian,
    approx_differentiation(
      function(t) as.numeric(concat_groups(t)),
      c(0.5, 2)
    ),
    tolerance = 1e-6
  )
})

# ---- a shaped tangent contributes one Jacobian column ------------------------
#
# `auto_differentiation()` assembles the Jacobian with `do.call(cbind, ...)`, so
# each column has to arrive as one entry per output. A scalar pair can carry a
# shaped tangent, which `theta[k] * X` produces and which the `dimnames<-`
# setter records labels on, and handing that back unflattened makes `cbind()`
# read a matrix as several columns: the Jacobian comes back with the payload's
# rows and one column per parameter for every payload column, and any labels the
# function recorded ride into it. The tangent-array branch flattens with
# `as.vector()` for exactly this reason, and the scalar-pair branch reads the
# same way.

test_that("extract_tangent_column flattens a shaped tangent on a scalar pair", {
  pair <- primal_tangent(
    base::matrix(c(1, 2, 3, 4), nrow = 2),
    base::matrix(c(5, 6, 7, 8), nrow = 2)
  )
  column <- extract_tangent_column(pair)
  expect_null(dim(column))
  expect_equal(column, c(5, 6, 7, 8))
})

test_that("a shaped-tangent result gives one Jacobian row per output", {
  X <- base::matrix(c(1, 2, 3, 4), nrow = 2)
  scale_design <- function(theta) theta[1] * X
  jacobian <- auto_differentiation(c(3, 5), scale_design)
  # Flattening the same result with c() before returning it reaches the
  # tangent-array branch, which already reports one column per parameter, so it
  # is the reference for the shape as well as for the values. The two are
  # compared axis by axis and then value by value, because comparing two numeric
  # matrices of different shapes directly reports the shapes through an error
  # from inside the comparison rather than as a difference.
  flattened <- auto_differentiation(c(3, 5), function(theta) c(theta[1] * X))
  expect_equal(dim(jacobian), c(4L, 2L))
  expect_equal(dim(jacobian), dim(flattened))
  expect_equal(as.vector(jacobian), as.vector(flattened))
  expect_equal(as.vector(jacobian), c(1, 2, 3, 4, 0, 0, 0, 0))
})

test_that("labels on a shaped tangent stay out of the Jacobian", {
  # A labeled psi return reports the unnamed bread the numeric pass reports,
  # which holds because a Jacobian column is read off the tangent slot through
  # as.vector() or off its row sums, and both drop the labels. The scalar-pair
  # branch has to drop them too, or a function that labels its rows reports a
  # labeled Jacobian while the same function labeling nothing does not.
  label_rows <- function(theta) {
    v <- theta[1] * base::matrix(c(1, 2), ncol = 1)
    base::`rownames<-`(v, c("a", "b"))
  }
  jacobian <- auto_differentiation(c(2, 7), label_rows)
  expect_equal(dim(jacobian), c(2L, 2L))
  expect_null(dimnames(jacobian))
  expect_equal(jacobian, base::cbind(c(1, 2), c(0, 0)))
})

test_that("names on a tangent stay out of the Jacobian", {
  # Multiplying a scalar pair through a named vector carries the names onto both
  # slots, which is the same leak without a matrix in sight.
  named_outputs <- function(theta) theta[1] * c(a = 1, b = 2)
  jacobian <- auto_differentiation(c(2.5, 3.75), named_outputs)
  expect_null(dimnames(jacobian))
  expect_equal(jacobian, base::cbind(c(1, 2), c(0, 0)))
})

test_that("the already flat tangent shapes report the same column", {
  expect_equal(extract_tangent_column(primal_tangent(2, 3)), 3)
  expect_equal(
    extract_tangent_column(primal_tangent(c(1, 2, 3), c(4, 5, 6))),
    c(4, 5, 6)
  )
  expect_equal(
    extract_tangent_column(primal_tangent_array(
      base::matrix(c(1, 2, 3, 4), nrow = 2),
      base::matrix(c(9, 8, 7, 6), nrow = 2)
    )),
    c(9, 8, 7, 6)
  )
  expect_equal(
    extract_tangent_column(
      primal_tangent_vector(list(primal_tangent(1, 1), primal_tangent(2, 0)))
    ),
    c(1, 0)
  )
})

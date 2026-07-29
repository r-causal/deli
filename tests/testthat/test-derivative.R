# ---- linear function ---------------------------------------------------------

test_that("approx_differentiation computes Jacobian for linear function", {
  # f(x) = c(2*x[1] + 3*x[2], x[1] - x[2])
  # Jacobian: [[2, 3], [1, -1]] (constant, exact for all methods)
  f <- function(x) c(2 * x[1] + 3 * x[2], x[1] - x[2])
  theta <- c(1, 1)
  expected <- matrix(c(2, 1, 3, -1), nrow = 2, ncol = 2)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(result, expected, tolerance = 1e-5)
  }
})

# ---- quadratic function -----------------------------------------------------

test_that("approx_differentiation computes Jacobian for quadratic function", {
  # f(x) = c(x[1]^2, x[2]^2)
  # At theta = c(3, 4): Jacobian = [[6, 0], [0, 8]]
  f <- function(x) c(x[1]^2, x[2]^2)
  theta <- c(3, 4)
  expected <- matrix(c(6, 0, 0, 8), nrow = 2, ncol = 2)

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result, expected, tolerance = 1e-5)
})

test_that("quadratic Jacobian has zero off-diagonal elements", {
  f <- function(x) c(x[1]^2, x[2]^2)
  theta <- c(3, 4)

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 2], 0, tolerance = 1e-5)
  expect_equal(result[2, 1], 0, tolerance = 1e-5)
})

# ---- scalar function ---------------------------------------------------------

test_that("approx_differentiation handles scalar function f(x) = x^3", {
  # f(x) = x^3, derivative = 3x^2 = 12 at x = 2
  f <- function(x) x^3
  theta <- 2

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 1], 12, tolerance = 1e-5)
})

# ---- method agreement --------------------------------------------------------

test_that("central, forward, and backward methods agree for smooth functions", {
  f <- function(x) c(sin(x[1]) + x[2]^2, x[1] * x[2])
  theta <- c(1.5, 2.5)

  result_central <- approx_differentiation(
    func = f,
    theta = theta,
    method = "capprox"
  )
  result_forward <- approx_differentiation(
    func = f,
    theta = theta,
    method = "fapprox"
  )
  result_backward <- approx_differentiation(
    func = f,
    theta = theta,
    method = "bapprox"
  )

  expect_equal(result_central, result_forward, tolerance = 1e-5)
  expect_equal(result_central, result_backward, tolerance = 1e-5)
})

# ---- dx parameter ------------------------------------------------------------

test_that("smaller dx gives more accurate results for smooth functions", {
  # f(x) = x^2, derivative at x = 5 is 10
  f <- function(x) x^2
  theta <- 5
  exact <- 10

  result_large <- approx_differentiation(
    func = f,
    theta = theta,
    method = "fapprox",
    dx = 1e-3
  )
  result_small <- approx_differentiation(
    func = f,
    theta = theta,
    method = "fapprox",
    dx = 1e-9
  )

  error_large <- abs(result_large[1, 1] - exact)
  error_small <- abs(result_small[1, 1] - exact)
  expect_true(error_small < error_large)
})

# ---- exponential function ----------------------------------------------------

test_that("approx_differentiation computes derivative of exp(x) at x = 1", {
  f <- function(x) exp(x)
  theta <- 1

  result <- approx_differentiation(func = f, theta = theta, method = "capprox")
  expect_equal(result[1, 1], exp(1), tolerance = 1e-5)
})

# ---- column-matrix-returning function ----------------------------------------

test_that("approx_differentiation handles a func that returns a column matrix", {
  # f(th) = X %*% th returns an n-by-1 column matrix in R (any %*% does).
  # The Jacobian of a linear map is X itself, an n-by-p matrix.
  X <- matrix(c(0.3, -1.1, 0.7, 2.0, -0.5, 1.4), nrow = 3, ncol = 2)
  f <- function(th) X %*% th
  theta <- c(0.5, -0.3)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(3L, 2L))
    expect_equal(unname(result), X, tolerance = 1e-5)
  }
})

test_that("approx_differentiation returns n-by-p for vector-valued func", {
  # A plain vector return of length n with p parameters yields an n-by-p matrix.
  f <- function(x) c(2 * x[1] + 3 * x[2], x[1] - x[2], x[1])
  theta <- c(1, 1)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(3L, 2L))
  }
})

test_that("approx_differentiation returns 1-by-1 for scalar func of one parameter", {
  f <- function(x) x^3
  theta <- 2

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(1L, 1L))
    expect_equal(result[1, 1], 12, tolerance = 1e-5)
  }
})

test_that("approx_differentiation handles a length-1 vector return", {
  # A length-1 vector output with two parameters yields a 1-by-2 gradient row.
  f <- function(x) x[1] * x[2]
  theta <- c(3, 5)

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_equal(dim(result), c(1L, 2L))
    expect_equal(result[1, 1], 5, tolerance = 1e-5)
    expect_equal(result[1, 2], 3, tolerance = 1e-5)
  }
})

# ---- all methods run without error -------------------------------------------

test_that("all three methods run without error", {
  f <- function(x) c(x[1] + x[2], x[1] * x[2])
  theta <- c(2, 3)

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_no_error(
      approx_differentiation(func = f, theta = theta, method = method)
    )
  }
})

# ---- step size against parameter magnitude -----------------------------------
#
# `dx` is an absolute perturbation, so at a large enough parameter magnitude the
# requested step falls below the spacing of the doubles surrounding `theta`.
# `theta + dx` then rounds back to `theta`, the difference quotient loses its
# significant digits, and past roughly 1.7e7 it collapses to exactly zero. The
# step is therefore floored at the resolution of each parameter, which leaves
# `dx` itself in force wherever it is representable. The two regimes are pinned
# separately because a fix that only rescues the collapsed one leaves a silent
# error of several percent standing in the degraded one.

# A well-scaled estimating function at parameter magnitude `s`: the value grows
# with the data while the derivative stays order one, which is the shape of a
# real estimating equation fitted to data on that scale.
scaled_nonlinear <- function(s) {
  function(x) s * log1p((x / s)^2)
}

# d/dx of the above, in scaled units u = x / s.
scaled_nonlinear_deriv <- function(u) 2 * u / (1 + u^2)

# The floor guarantees the step spans enough of the parameter's resolution to
# reproduce four significant digits, so the surviving error stays a few parts in
# a hundred thousand at any magnitude. The measured worst case over the eighteen
# combinations of function, magnitude, and method that the two tests below run
# is 8.55e-05; the tolerance leaves room for platform variation.
step_floor_tolerance <- 5e-3

test_that("a dx near the resolution of theta gives the right derivative", {
  # At these magnitudes the requested step still perturbs `theta`, but it spans
  # only 550, 69, and 1.07 representable values, so the quotient silently loses
  # most of its digits. Without the floor the error over these nine combinations
  # runs from 2.59e-04 to 1.
  u <- 1.3
  expected <- scaled_nonlinear_deriv(u)

  for (method in c("capprox", "fapprox", "bapprox")) {
    for (s in c(1e4, 1e5, 5e6)) {
      result <- approx_differentiation(
        func = scaled_nonlinear(s),
        theta = s * u,
        method = method,
        dx = 1e-9
      )
      expect_equal(
        result[1, 1],
        expected,
        tolerance = step_floor_tolerance,
        label = paste0(method, " at |theta| = ", s * u)
      )
    }
  }
})

test_that("a dx below the resolution of theta gives the right derivative", {
  # Past roughly 1.7e7 the requested step is lost entirely: `theta + dx` is
  # `theta`, every difference is zero, and without the floor the derivative is
  # exactly zero with no warning of any kind.
  u <- 1.3
  expected <- scaled_nonlinear_deriv(u)

  for (method in c("capprox", "fapprox", "bapprox")) {
    for (s in c(5e8, 5e9, 5e10)) {
      result <- approx_differentiation(
        func = scaled_nonlinear(s),
        theta = s * u,
        method = method,
        dx = 1e-9
      )
      expect_false(result[1, 1] == 0)
      expect_equal(
        result[1, 1],
        expected,
        tolerance = step_floor_tolerance,
        label = paste0(method, " at |theta| = ", s * u)
      )
    }
  }
})

test_that("each Jacobian column is divided by its own parameter's step", {
  # Once the steps differ from one another, every difference has to meet the
  # step of the parameter that produced it. Here the two parameters are nine
  # orders of magnitude apart, so the first keeps `dx` itself while the floor
  # raises the second to 6.66e-3: an entry divided by the other parameter's
  # step is wrong by that ratio, 6.66e6. The function returns three values for
  # two parameters so that a transposed Jacobian cannot satisfy the assertion
  # either.
  f <- function(th) {
    a <- th[1]
    u <- th[2] / 1e9
    c(a^3 * u, sin(a) + u^2, exp(a / 4) * log(u))
  }
  theta <- c(2, 3e9)
  a <- theta[1]
  u <- theta[2] / 1e9
  expected <- matrix(
    c(
      3 * a^2 * u,
      cos(a),
      exp(a / 4) * log(u) / 4,
      a^3 / 1e9,
      2 * u / 1e9,
      exp(a / 4) / theta[2]
    ),
    nrow = 3,
    ncol = 2
  )

  for (method in c("capprox", "fapprox", "bapprox")) {
    result <- approx_differentiation(func = f, theta = theta, method = method)
    expect_identical(dim(result), c(3L, 2L))
    # Relative, entry by entry: the two columns of the Jacobian are themselves
    # nine orders of magnitude apart, so an absolute comparison would leave the
    # second column unconstrained. The measured worst case is 1.29e-05.
    expect_lt(
      max(abs(result / expected - 1)),
      1e-3,
      label = paste0("worst relative error, ", method)
    )
  }
})

test_that("the perturbation is exactly dx where dx is representable at theta", {
  # The floor may not disturb a well-scaled fit. These magnitudes are the ones
  # the package's own reference fixtures reach, so the quotient must be the one
  # a step of exactly `dx` produces, to the last bit.
  f <- function(x) exp(x / 1e3) + x^2 / 1e4
  dx <- 1e-9

  for (theta in c(2, 37.23, 183.6)) {
    by_hand <- (f(theta + dx) - f(theta - dx)) / (2 * dx)
    result <- approx_differentiation(
      func = f,
      theta = theta,
      method = "capprox",
      dx = dx
    )
    expect_identical(result[1, 1], by_hand)
  }
})

test_that("approx_differentiation honors a deliberately large dx", {
  # `ee_percentile()` is differentiated with `dx = 1` to smooth an indicator,
  # matching the recipe Python delicatessen documents, so a large `dx` has to
  # stay the absolute step the caller asked for rather than being rescaled.
  f <- function(x) x^2
  theta <- 5.4

  result <- approx_differentiation(
    func = f,
    theta = theta,
    method = "fapprox",
    dx = 1
  )
  expect_identical(result[1, 1], (f(theta + 1) - f(theta)) / 1)
  # A step that large is deliberately nowhere near the true derivative of 10.8.
  expect_equal(result[1, 1], 11.8)
})

test_that("a genuinely flat function differentiates to zero without warning", {
  # A step function has a finite-difference derivative of exactly zero because
  # the function does not change, not because the step was lost. The two cases
  # must not be conflated: this one is correct and must stay quiet.
  f <- function(x) as.numeric(x > 0)
  theta <- 2.5

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_no_warning({
      result <- approx_differentiation(
        func = f,
        theta = theta,
        method = method
      )
    })
    expect_identical(result[1, 1], 0)
  }
})

# ---- a difference lost against the function values ---------------------------
#
# The floor above rescues a step lost against a large `theta`. The other way a
# difference vanishes is untouched by it. Where the step is applied exactly but
# the function values are large, the change it produces falls below the spacing
# of the doubles holding them: both evaluations round to the same value and the
# quotient is exactly zero. Nothing about the step is wrong, so no reading of the
# step can see this, and the result is indistinguishable from the genuinely flat
# function pinned above. What separates the two is the significance of the
# difference against the magnitude of the values it was taken between. A
# difference carrying none of it, at a step fine enough that a real derivative
# would have shown, has to be reported rather than returned as a zero.

test_that("a difference lost against large function values is reported", {
  # The values are about 1e11, whose neighboring doubles are 1.53e-05 away,
  # while the default step moves the function by 4e-07. Both evaluations return
  # the same double, so the derivative comes back as zero rather than 200.
  f <- function(x) 1e11 + 200 * x

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_warning(
      approx_differentiation(func = f, theta = 0, method = method, dx = 1e-9),
      class = "deli_finite_difference_lost",
      label = method
    )
  }
})

test_that("the same values differentiate quietly at a step that resolves them", {
  # Nothing about the function has changed, only the step: 1e-03 moves it by
  # 0.4, some 26,000 times the spacing of the values, so the difference keeps
  # its digits and the derivative is the right one.
  f <- function(x) 1e11 + 200 * x

  for (method in c("capprox", "fapprox", "bapprox")) {
    expect_no_warning({
      result <- approx_differentiation(
        func = f,
        theta = 0,
        method = method,
        dx = 1e-3
      )
    })
    expect_equal(result[1, 1], 200, tolerance = 1e-4, label = method)
  }
})

# Exact-mode (forward-mode autodiff) sandwich variance for the pooled logistic
# survival estimating equation. Under deriv_method = "exact", the covariate
# log-odds carry tangents, and ee_plogit tiles them across the unit-time
# intervals with a replication and sums the interval residuals down the columns.
# Those reshaping operations must preserve the tangent, so these tests pin two
# invariants: the exact-mode bread, variance, and intervals match the Python
# Delicatessen reference computed with deriv_method='exact', and they agree with
# the finite-difference derivative on this smooth problem.
#
# Tolerance note: the internal exact-versus-capprox checks use 1e-6 because a
# central difference with dx = 1e-9 carries a rounding floor near 1e-7. The
# fixture comparison also uses 1e-6: the reference theta comes from Python's
# Levenberg-Marquardt solver while the R fit uses rootSolve, and the small
# difference in the solved root propagates into the bread at that scale.

# ---- fixture parity: exact-mode results match Python ------------------------

test_that("ee_plogit exact mode matches the Python reference", {
  ref <- load_fixture("ee_plogit_exact")
  psi <- function(theta) {
    ee_plogit(theta, X = ref$X, time = ref$time, event = ref$event)
  }
  m <- MEstimator(stacked_equations = psi, init = ref$init)
  m <- estimate(m, deriv_method = "exact")
  expect_python_match(m, "ee_plogit_exact", tolerance = 1e-6)
})

# ---- internal consistency: exact agrees with the finite difference ----------

# The exact bread must also agree with the central-difference bread on the same
# fit, which isolates the tangent-aware replication and column summation from the
# solver. The point estimate is identical between the two fits (the derivative
# method only enters the bread), so this isolates the derivative computation.
test_that("exact and capprox derivatives agree for ee_plogit", {
  ref <- load_fixture("ee_plogit_exact")
  psi <- function(theta) {
    ee_plogit(theta, X = ref$X, time = ref$time, event = ref$event)
  }
  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  m_cap <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "capprox"
  )
  expect_equal(unname(m_exact@bread), unname(m_cap@bread), tolerance = 1e-6)
  expect_equal(
    sqrt(diag(m_exact@variance)),
    sqrt(diag(m_cap@variance)),
    tolerance = 1e-6
  )
})

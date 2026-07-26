# Exact-mode (forward-mode autodiff) sandwich variance for the causal and
# measurement estimating equations. Under deriv_method = "exact", theta enters an
# estimating function as a tangent-carrying object, and every linear predictor
# X %*% theta must flatten to a vector without dropping the tangent. These tests
# pin two invariants: the exact-mode theta, bread, variance, and intervals match
# the Python Delicatessen reference computed with deriv_method='exact', and they
# agree with the finite-difference derivative on these problems.
#
# Tolerance note: the internal exact-versus-capprox checks use 1e-6 because a
# central difference with dx = 1e-9 carries a rounding floor near 1e-7, so exact
# and capprox derivatives cannot be expected to agree more tightly than that. The
# fixture comparisons also use 1e-6: the reference theta comes from Python's
# Levenberg-Marquardt solver while the R fit uses rootSolve, and the small
# difference in the solved root propagates into the bread at that scale.

# Compare the exact bread and standard errors against the central-difference
# bread and standard errors on a single fit. The point estimate is identical
# between the two fits (the derivative method only enters the bread), so this
# isolates the derivative computation. The standard-error magnitudes here are all
# above 0.03, so the 1e-6 comparison is not vacuous.
expect_exact_matches_capprox_cm <- function(psi, init) {
  m_exact <- estimate(
    MEstimator(stacked_equations = psi, init = init),
    deriv_method = "exact"
  )
  m_cap <- estimate(
    MEstimator(stacked_equations = psi, init = init),
    deriv_method = "capprox"
  )
  expect_equal(unname(m_exact@bread), unname(m_cap@bread), tolerance = 1e-6)
  expect_equal(
    sqrt(diag(m_exact@variance)),
    sqrt(diag(m_cap@variance)),
    tolerance = 1e-6
  )
}

# ---- fixture parity: exact-mode results match Python ------------------------

test_that("ee_ipw exact mode matches the Python reference", {
  ref <- load_fixture("ee_ipw_exact")
  psi <- function(theta) ee_ipw(theta, y = ref$y, A = ref$a, W = ref$W)
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_ipw_exact", tolerance = 1e-6)
})

test_that("ee_gformula exact mode matches the Python reference", {
  ref <- load_fixture("ee_gformula_exact")
  psi <- function(theta) {
    ee_gformula(theta, y = ref$y, X = ref$X, X1 = ref$X1, X0 = ref$X0)
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_gformula_exact", tolerance = 1e-6)
})

test_that("ee_aipw exact mode matches the Python reference", {
  ref <- load_fixture("ee_aipw_exact")
  psi <- function(theta) {
    ee_aipw(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      X = ref$X,
      X1 = ref$X1,
      X0 = ref$X0
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_aipw_exact", tolerance = 1e-6)
})

test_that("ee_gestimation_snmm exact mode matches the Python reference", {
  ref <- load_fixture("ee_gestimation_exact")
  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      model = ref$model
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_gestimation_exact", tolerance = 1e-6)
})

# A single-column structural mean model (V is the intercept alone) makes the SMM
# parameter subset theta[1] a lone scalar pair rather than a length-1 tangent
# vector, so the matrix product (V*A) %*% phi carries tangents through the scalar
# branch of the product. This pins the exact-mode results for that branch.
test_that("ee_gestimation_snmm single-column V exact mode matches the Python reference", {
  ref <- load_fixture("ee_gestimation_single_v_exact")
  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      model = ref$model
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_gestimation_single_v_exact", tolerance = 1e-6)
})

# The efficient g-estimator fits an outcome model E[h(phi) | X] where the
# response h(phi) = Y - (V*A) phi is theta-derived, so it carries tangents under
# exact mode and reaches ee_regression as a tangent-carrying y argument. This
# pins the exact-mode results for that routing.
test_that("ee_gestimation_snmm efficient exact mode matches the Python reference", {
  ref <- load_fixture("ee_gestimation_efficient_exact")
  psi <- function(theta) {
    ee_gestimation_snmm(
      theta,
      y = ref$y,
      A = ref$a,
      W = ref$W,
      V = ref$V,
      X = ref$X,
      model = ref$model
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_gestimation_efficient_exact", tolerance = 1e-6)
})

# Two-stage least squares stacks predicted A (theta-derived, tangent-carrying)
# into the second-stage design and routes it through ee_regression as the X
# argument. The with-W case gives a multi-column design; the no-W case gives a
# single tangent-carrying column with a scalar second-stage coefficient.
test_that("ee_2sls with W exact mode matches the Python reference", {
  ref <- load_fixture("ee_2sls_exact_w")
  psi <- function(theta) {
    ee_2sls(theta, y = ref$y, A = ref$a, Z = ref$Z, W = ref$W)
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_2sls_exact_w", tolerance = 1e-6)
})

test_that("ee_2sls without W exact mode matches the Python reference", {
  ref <- load_fixture("ee_2sls_exact_now")
  psi <- function(theta) {
    ee_2sls(theta, y = ref$y, A = ref$a, Z = ref$Z)
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_2sls_exact_now", tolerance = 1e-6)
})

test_that("ee_rogan_gladen_extended exact mode matches the Python reference", {
  ref <- load_fixture("ee_rogan_gladen_extended_exact")
  psi <- function(theta) {
    ee_rogan_gladen_extended(
      theta,
      y = ref$y,
      y_star = ref$y_star,
      r = ref$r,
      X = ref$X
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_rogan_gladen_extended_exact", tolerance = 1e-6)
})

test_that("ee_regression_calibration exact mode matches the Python reference", {
  ref <- load_fixture("ee_regression_calibration_exact")
  psi <- function(theta) {
    ee_regression_calibration(
      theta,
      beta = ref$beta,
      a = ref$a,
      a_star = ref$a_star,
      r = ref$r
    )
  }
  m <- estimate(
    MEstimator(stacked_equations = psi, init = ref$init),
    deriv_method = "exact"
  )
  expect_python_match(m, "ee_regression_calibration_exact", tolerance = 1e-6)
})

# ---- internal consistency: exact agrees with the finite difference ----------

test_that("exact and capprox derivatives agree for the causal equations", {
  ipw <- load_fixture("ee_ipw_exact")
  expect_exact_matches_capprox_cm(
    function(t) ee_ipw(t, y = ipw$y, A = ipw$a, W = ipw$W),
    ipw$init
  )
  gf <- load_fixture("ee_gformula_exact")
  expect_exact_matches_capprox_cm(
    function(t) ee_gformula(t, y = gf$y, X = gf$X, X1 = gf$X1, X0 = gf$X0),
    gf$init
  )
  aipw <- load_fixture("ee_aipw_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_aipw(
        t,
        y = aipw$y,
        A = aipw$a,
        W = aipw$W,
        X = aipw$X,
        X1 = aipw$X1,
        X0 = aipw$X0
      )
    },
    aipw$init
  )
  ge <- load_fixture("ee_gestimation_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_gestimation_snmm(
        t,
        y = ge$y,
        A = ge$a,
        W = ge$W,
        V = ge$V,
        model = ge$model
      )
    },
    ge$init
  )
  gv <- load_fixture("ee_gestimation_single_v_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_gestimation_snmm(
        t,
        y = gv$y,
        A = gv$a,
        W = gv$W,
        V = gv$V,
        model = gv$model
      )
    },
    gv$init
  )
  ge_eff <- load_fixture("ee_gestimation_efficient_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_gestimation_snmm(
        t,
        y = ge_eff$y,
        A = ge_eff$a,
        W = ge_eff$W,
        V = ge_eff$V,
        X = ge_eff$X,
        model = ge_eff$model
      )
    },
    ge_eff$init
  )
  tsw <- load_fixture("ee_2sls_exact_w")
  expect_exact_matches_capprox_cm(
    function(t) ee_2sls(t, y = tsw$y, A = tsw$a, Z = tsw$Z, W = tsw$W),
    tsw$init
  )
  tsn <- load_fixture("ee_2sls_exact_now")
  expect_exact_matches_capprox_cm(
    function(t) ee_2sls(t, y = tsn$y, A = tsn$a, Z = tsn$Z),
    tsn$init
  )
})

test_that("exact and capprox derivatives agree for the measurement equations", {
  rg <- load_fixture("ee_rogan_gladen_extended_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_rogan_gladen_extended(
        t,
        y = rg$y,
        y_star = rg$y_star,
        r = rg$r,
        X = rg$X
      )
    },
    rg$init
  )
  rc <- load_fixture("ee_regression_calibration_exact")
  expect_exact_matches_capprox_cm(
    function(t) {
      ee_regression_calibration(
        t,
        beta = rc$beta,
        a = rc$a,
        a_star = rc$a_star,
        r = rc$r
      )
    },
    rc$init
  )
})

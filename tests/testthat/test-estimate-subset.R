# ---- Subset solving only solves specified parameters -------------------------

test_that("subset solving only solves for specified parameters", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, only solve for theta[2] (variance)
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # theta[2] should be solved to the reference variance

  expect_equal(unname(m@theta[2]), ref$theta[2], tolerance = 1e-4)
})

# ---- Non-subset parameters remain at init values ----------------------------

test_that("non-subset parameters remain at init values", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, only solve for theta[2]
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # theta[1] should remain at the init value (the true mean)
  expect_equal(unname(m@theta[1]), true_mean)
})

# ---- Results match when subset includes all parameters -----------------------

test_that("results match when subset includes all parameters", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Estimate without subset
  m_full <- MEstimator(stacked_equations = psi, init = c(0, 1))
  m_full <- estimate(m_full)

  # Estimate with subset = all parameters
  m_subset <- MEstimator(
    stacked_equations = psi,
    init = c(0, 1),
    subset = c(1L, 2L)
  )
  m_subset <- estimate(m_subset)

  expect_equal(m_subset@theta, m_full@theta, tolerance = 1e-6)
})

# ---- Subset with mean+variance EE -------------------------------------------

test_that("subset with mean+variance EE matches reference values", {
  ref <- load_fixture("ee_mean_variance")
  y <- ref$y

  psi <- function(theta) {
    rbind(
      y - theta[1],
      (y - theta[1])^2 - theta[2]
    )
  }

  # Fix theta[1] at the true mean, solve only for variance (theta[2])
  true_mean <- ref$theta[1]
  m <- MEstimator(
    stacked_equations = psi,
    init = c(true_mean, 1),
    subset = 2L
  )
  m <- estimate(m)

  # The full theta should combine the fixed mean and the solved variance
  expect_equal(unname(m@theta[1]), true_mean)
  expect_equal(unname(m@theta[2]), ref$theta[2], tolerance = 1e-4)

  # bread, meat, and variance should still be computed
  expect_true(!is.null(m@bread))
  expect_true(!is.null(m@meat))
  expect_true(!is.null(m@variance))
})

# ---- The fitted GMM path under a subset --------------------------------------
#
# The wrapper tests pin that `subset` reaches the estimator. What the fitted path
# does with it is separate, and three things about it are documented rather than
# obvious.
#
# The variance estimator ignores `subset` (see gmm_estimate()'s `@param
# subset`): the bread and the meat are built over every parameter and every
# moment condition, so the covariance is the full one at the point the fit
# returned, whatever the minimizer was allowed to move.
#
# The objective is a quadratic form in every moment condition, so the conditions
# outside the subset are still summed in and still pull on the free parameters.
# A subset GMM fit is therefore not the fit of the subset equations on their
# own, which is what MEstimator() returns from the same stack.
#
# And a subset fit is not judged on its moments. The parameters outside the
# subset are held at their initial values rather than estimated, which neither
# the vanishing-moments reading of a just-identified fit nor the reference
# distribution of the J-statistic allows for, so both are left unmade.

subset_gmm_data <- function() {
  set.seed(31)
  n <- 120
  x <- stats::rnorm(n)
  list(x = x, y = 1 + 2 * x + stats::rnorm(n), n = n)
}

subset_gmm_just <- function(d) {
  function(theta) {
    r <- d$y - theta[1] - theta[2] * d$x
    rbind(r, r * d$x)
  }
}

subset_gmm_over <- function(d) {
  function(theta) {
    r <- d$y - theta[1] - theta[2] * d$x
    rbind(r, r * d$x, r * d$x^2)
  }
}

test_that("a subset GMM fit reports the variance of the whole parameter vector", {
  d <- subset_gmm_data()
  psi <- subset_gmm_just(d)
  g <- gmm_estimate(psi, init = c(0, 1), subset = 1L)

  # The frozen parameter keeps its starting value and still gets a row and a
  # column of the covariance.
  expect_equal(unname(coef(g))[[2L]], 1)
  expect_equal(dim(vcov(g)), c(2L, 2L))
  expect_equal(dim(g@bread), c(2L, 2L))
  expect_equal(dim(g@meat), c(2L, 2L))

  # The variance estimator ignores `subset`, so it is the sandwich of the whole
  # stack at the point the fit returned and nothing about it is restricted.
  expect_equal(
    unname(vcov(g)),
    compute_sandwich(psi, theta = unname(coef(g))) / d$n,
    tolerance = 1e-8
  )
})

test_that("a subset GMM fit is not the fit of the subset equations alone", {
  d <- subset_gmm_data()
  psi <- subset_gmm_just(d)

  gmm <- gmm_estimate(psi, init = c(0, 1), subset = 1L)
  m <- m_estimate(psi, init = c(0, 1), subset = 1L)

  # Both hold the slope where it started, and they disagree about the intercept:
  # the minimizer weighs the second moment condition in even though nothing it
  # may move appears only there, while the root-finder sets the second equation
  # aside with the parameter it estimates.
  expect_equal(unname(coef(gmm))[[2L]], 1)
  expect_equal(unname(coef(m))[[2L]], 1)
  expect_false(isTRUE(all.equal(
    unname(coef(gmm))[[1L]],
    unname(coef(m))[[1L]]
  )))
})

test_that("a just-identified subset fit leaves the weight matrix at the identity", {
  d <- subset_gmm_data()
  g <- gmm_estimate(subset_gmm_just(d), init = c(0, 1), subset = 1L)

  # The two-step update runs only for an over-identified system, so nothing
  # moves the weight matrix here, and it is indexed by the moment conditions
  # rather than by the parameters, so it carries no labels.
  expect_equal(g@weight_matrix, diag(2))
  expect_null(dimnames(g@weight_matrix))
})

test_that("an over-identified subset fit still updates the weight matrix", {
  skip_if_not_installed("MASS")
  d <- subset_gmm_data()
  g <- gmm_estimate(
    subset_gmm_over(d),
    init = c(0, 2),
    subset = 1L,
    overid_maxiter = 200L
  )

  # Three moment conditions for two parameters, so the update runs even though
  # only one parameter is free to move within the objective.
  expect_equal(dim(g@weight_matrix), c(3L, 3L))
  expect_false(isTRUE(all.equal(g@weight_matrix, diag(3))))
  expect_null(dimnames(g@weight_matrix))

  # The bread is rectangular, indexed by the moment conditions down its rows and
  # by every parameter across its columns, and the covariance it pseudo-inverts
  # into is square in the parameters.
  expect_equal(dim(g@bread), c(3L, 2L))
  expect_equal(dim(g@meat), c(3L, 3L))
  expect_equal(dim(vcov(g)), c(2L, 2L))
  expect_equal(
    unname(vcov(g)),
    compute_sandwich(subset_gmm_over(d), theta = unname(coef(g))) / d$n,
    tolerance = 1e-8
  )
})

test_that("a subset fit records no J-statistic even when over-identified", {
  skip_if_not_installed("MASS")
  d <- subset_gmm_data()
  g <- gmm_estimate(
    subset_gmm_over(d),
    init = c(0, 2),
    subset = 1L,
    overid_maxiter = 200L
  )

  expect_null(g@j_statistic)
  expect_identical(generics::glance(g)$j_statistic, NA_real_)
  expect_identical(generics::glance(g)$j_df, NA_integer_)
})

test_that("a subset fit leaves its moments unjudged", {
  d <- subset_gmm_data()
  # The slope is frozen far from where the data put it, so the moment
  # conditions cannot vanish at the point the fit returns. A fit of the same
  # stack without a subset would be judged on exactly that and would report it.
  psi <- subset_gmm_just(d)

  expect_no_warning(gmm_estimate(psi, init = c(0, -5), subset = 1L))

  frozen <- gmm_estimate(psi, init = c(0, -5), subset = 1L)
  moments <- rowSums(psi(unname(coef(frozen))))
  expect_gt(max(abs(moments)), 1)
})

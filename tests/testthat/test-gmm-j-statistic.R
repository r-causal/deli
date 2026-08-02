# Hansen's J-statistic for over-identified GMM fits.
#
# A just-identified fit has as many moment conditions as parameters, so the
# moments vanish at a solution and the size of what is left over is a reading of
# whether the fit succeeded. estimate() takes that reading and warns on it. An
# over-identified fit has no such reading available, because no value of the
# parameters drives every moment to zero and a residual moment is expected
# rather than diagnostic. That left over-identified fits with no moment-quality
# instrument at all: a badly specified one could sit at a summed moment in the
# thousands and say nothing about it.
#
# J = n * gbar' W gbar at the optimum is the reading that is available. Under
# correct specification it is asymptotically chi-squared on m - p degrees of
# freedom, m moment conditions against p parameters, so its size can be judged
# against a reference distribution rather than against the scale of the data.
#
# The weight matrix is what makes J comparable across problems, so J means what
# it is supposed to mean only once the two-step update has settled on the
# efficient weight matrix. At the identity weight matrix, which is where
# `overid_maxiter = 0` leaves a fit, J is an unstandardized sum of squared
# moments and carries no reference distribution at all. Every fit below
# therefore states its own budget rather than relying on the default, so these
# tests say nothing about what the default is and do not move when it does.

settled_budget <- 200L

j_reference <- function(fit, psi) {
  moments <- rowSums(psi(fit@theta)) / fit@n_obs
  as.numeric(
    fit@n_obs * crossprod(moments, fit@weight_matrix %*% moments)
  )
}

j_df <- function(fit) {
  nrow(fit@weight_matrix) - length(fit@theta)
}

j_p_value <- function(fit) {
  stats::pchisq(fit@j_statistic, j_df(fit), lower.tail = FALSE)
}

# ---- well-specified over-identified fits --------------------------------------
#
# Every fit in this block is correctly specified, so its J is a draw from the
# reference distribution and its P-value is an ordinary one. Measured across the
# over-identified fits the test suite already holds, the P-values run from 0.16
# to 0.86; across two hundred seeds of a three-instrument IV sweep the smallest
# is 1.1e-3. Nothing correctly specified came anywhere near the threshold below.

# A Poisson count is identified twice over, by its mean and by its variance,
# which are equal. This is the system the GMMEstimator() examples fit.
make_poisson_psi <- function(seed = 42, n = 200) {
  set.seed(seed)
  counts <- stats::rpois(n, lambda = 3)
  function(theta) {
    rbind(counts - theta[1], (counts - theta[1])^2 - theta[1])
  }
}

# A normal mean with the variance condition pinned at the value that generated
# the data, so both conditions are true.
make_mean_variance_psi <- function(variance = 2, seed = 5, n = 200) {
  set.seed(seed)
  y <- stats::rnorm(n, mean = 3, sd = sqrt(2))
  function(theta) {
    rbind(y - theta[1], (y - theta[1])^2 - variance)
  }
}

# A linear model with a fourth moment condition from a variable the outcome does
# not depend on, so the extra condition is valid and the model is correct.
make_valid_extra_moment_psi <- function(seed = 13, n = 400) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  spare <- stats::rnorm(n)
  y <- 1 + 2 * x1 + 3 * x2 + stats::rnorm(n, sd = 0.5)
  design <- cbind(1, x1, x2)
  function(theta) {
    residual <- as.numeric(y - design %*% theta)
    rbind(residual, x1 * residual, x2 * residual, spare * residual)
  }
}

well_specified_overid_cases <- function() {
  list(
    poisson = list(psi = make_poisson_psi(), init = 1),
    mean_variance = list(psi = make_mean_variance_psi(), init = 2),
    valid_extra = list(psi = make_valid_extra_moment_psi(), init = c(0, 0, 0))
  )
}

test_that("an over-identified fit carries a J-statistic", {
  fit <- gmm_estimate(
    stacked_equations = make_poisson_psi(),
    init = 1,
    overid_maxiter = settled_budget
  )
  expect_true(is.numeric(fit@j_statistic))
  expect_length(fit@j_statistic, 1L)
  expect_true(is.finite(fit@j_statistic))
  expect_gt(fit@j_statistic, 0)
})

test_that("the J-statistic equals n times the objective at the optimum", {
  # An independent computation of n * gbar' W gbar from the fit's own parts.
  for (case in well_specified_overid_cases()) {
    fit <- gmm_estimate(
      stacked_equations = case$psi,
      init = case$init,
      overid_maxiter = settled_budget
    )
    expect_equal(fit@j_statistic, j_reference(fit, case$psi), tolerance = 1e-10)
  }
})

test_that("the J-statistic is a plain sum of squared moments at the identity weight matrix", {
  # `overid_maxiter = 0` leaves the identity in place, where the general formula
  # collapses to n * sum(gbar^2) and can be written out by hand.
  psi <- make_mean_variance_psi()
  fit <- gmm_estimate(stacked_equations = psi, init = 2, overid_maxiter = 0L)
  expect_identical(fit@weight_matrix, diag(2))
  moments <- rowSums(psi(fit@theta)) / fit@n_obs
  expect_equal(fit@j_statistic, fit@n_obs * sum(moments^2), tolerance = 1e-10)
})

test_that("a well-specified over-identified fit reports an unremarkable J", {
  # Measured P-values: 0.744 for the Poisson pair, 0.864 for the mean and
  # variance, 0.535 for the linear model with a valid extra moment.
  for (case in well_specified_overid_cases()) {
    fit <- expect_no_warning(gmm_estimate(
      stacked_equations = case$psi,
      init = case$init,
      overid_maxiter = settled_budget
    ))
    expect_identical(j_df(fit), 1L)
    expect_gt(j_p_value(fit), 0.05)
  }
})

# ---- just-identified fits have no J ------------------------------------------

test_that("a just-identified GMM fit has no J-statistic", {
  # m - p is zero, so there is no reference distribution and nothing to report.
  design <- cbind(1, mtcars$wt, mtcars$hp)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
  fit <- gmm_estimate(stacked_equations = psi, init = c(0, 0, 0))
  expect_identical(j_df(fit), 0L)
  expect_null(fit@j_statistic)
})

test_that("a one-equation just-identified GMM fit has no J-statistic", {
  y <- c(1, 2, 3, 4, 5)
  psi <- function(theta) matrix(y - theta[1], nrow = 1)
  fit <- gmm_estimate(stacked_equations = psi, init = 0)
  expect_null(fit@j_statistic)
})

test_that("an unfitted GMMEstimator has no J-statistic", {
  psi <- make_poisson_psi()
  expect_null(GMMEstimator(stacked_equations = psi, init = 1)@j_statistic)
})

test_that("an MEstimator has no J-statistic property", {
  # J belongs to the over-identified case, which M-estimation cannot reach.
  design <- cbind(1, mtcars$wt)
  psi <- function(theta) {
    ee_regression(theta, X = design, y = mtcars$mpg, model = "linear")
  }
  m <- m_estimate(stacked_equations = psi, init = c(0, 0))
  expect_false("j_statistic" %in% names(S7::props(m)))
})

# ---- a subset fit is not judged ----------------------------------------------
#
# `subset` holds the parameters outside it at their initial values rather than
# estimating them. The degrees of freedom the reference distribution counts are
# the moment conditions beyond the parameters, and a fit that estimated fewer
# parameters than it has is not the system those degrees of freedom describe, so
# no J is reported at all rather than one measured against the wrong
# distribution. This is the same reason the summed-moment reading estimate()
# takes for a just-identified fit passes over a subset fit.

test_that("an over-identified subset GMM fit has no J-statistic", {
  # Three moment conditions against two parameters, one of them held fixed at the
  # value that generated the data. The system is over-identified and correctly
  # specified, so nothing here would warn on its own.
  set.seed(29)
  n <- 300
  x1 <- stats::rnorm(n)
  spare <- stats::rnorm(n)
  y <- 1 + 2 * x1 + stats::rnorm(n, sd = 0.5)
  design <- cbind(1, x1)
  psi <- function(theta) {
    residual <- as.numeric(y - design %*% theta)
    rbind(residual, x1 * residual, spare * residual)
  }
  fit <- expect_no_warning(gmm_estimate(
    stacked_equations = psi,
    init = c(intercept = 0, slope = 2),
    subset = 1L,
    overid_maxiter = settled_budget
  ))
  # The count that would have been the degrees of freedom, so a system that
  # stopped being over-identified could not pass this quietly.
  expect_identical(j_df(fit), 1L)
  expect_null(fit@j_statistic)
  # And nothing in the printed summary either, which reads the property.
  summary_lines <- cli::cli_fmt(print(summary(fit)))
  expect_false(any(grepl("J-statistic", summary_lines, fixed = TRUE)))
})

# ---- the badly specified over-identified fit that used to stay quiet ---------
#
# The case this check exists for. Each fit below is an over-identified system
# whose moment conditions cannot all be true, and every one of them fits without
# a word: the summed-moment reading estimate() takes for a just-identified fit
# is not available here, so nothing looked at the moments at all.
#
# Measured J P-values, against the 0.16 to 0.86 the well-specified fits above
# report: 5.5e-31 for the omitted covariate, 1.0e-26 for the variance condition
# pinned at two and a half times the truth, and 9.6e-17 for the mean and
# variance of a sample whose spread the condition understates by a factor of
# thirty. The separation between the two groups is fourteen orders of magnitude
# at its narrowest, so the threshold has room to sit far from both.

# A linear model with a covariate left out of the design but kept as a moment
# condition, so that condition is false at every value of the parameters.
make_omitted_covariate_psi <- function(seed = 13, n = 400) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- 1 + 2 * x1 + 3 * x2 + stats::rnorm(n, sd = 0.5)
  design <- cbind(1, x1)
  function(theta) {
    residual <- as.numeric(y - design %*% theta)
    rbind(residual, x1 * residual, x2 * residual)
  }
}

# The mean and variance of a sample with a standard deviation of eight, against
# a variance condition pinned at two. This reproduces the reported case: the
# largest summed moment is 11914 and the fit says nothing about it.
make_understated_variance_psi <- function(seed = 21, n = 200) {
  set.seed(seed)
  y <- stats::rnorm(n, mean = 50, sd = 8)
  function(theta) {
    rbind(y - theta[1], (y - theta[1])^2 - 2)
  }
}

# The two fits below state a tolerance as well as a budget. This stack's weight
# update approaches its fixed point slowly, and how far into that approach a
# platform gets is where macOS and Linux part company: at 1e-6 the update
# settles in 97 passes on both, at 1e-7 in 151, and at the default of 1e-9
# within this budget on macOS but not on Linux. An update that has not settled
# leaves the fit with no J reading to judge, and an over-identification report
# in place of the one these tests are about, so the tolerance is stated rather
# than inherited. The statistic does not turn on the choice: J is 133.9852 at
# 1e-6 against 133.9848 at the default, on one degree of freedom either way.

test_that("a badly specified over-identified fit warns on its J-statistic", {
  cnd <- rlang::catch_cnd(
    gmm_estimate(
      stacked_equations = make_omitted_covariate_psi(),
      init = c(0, 0),
      overid_maxiter = settled_budget,
      overid_tolerance = 1e-6
    ),
    classes = "warning"
  )
  expect_s3_class(cnd, "rlang_warning")
  expect_s3_class(cnd, "deli_gmm_moments_rejected")
  expect_match(conditionMessage(cnd), "J-statistic")
})

test_that("the J warning reports the statistic it is judging", {
  # A warning that names no number leaves the reader unable to tell a marginal
  # failure from this one, which is off the scale of the reference distribution.
  cnd <- rlang::catch_cnd(
    gmm_estimate(
      stacked_equations = make_omitted_covariate_psi(),
      init = c(0, 0),
      overid_maxiter = settled_budget,
      overid_tolerance = 1e-6
    ),
    classes = "deli_gmm_moments_rejected"
  )
  # The measured statistic is 134.0 on one degree of freedom.
  expect_match(gsub("\\s+", " ", conditionMessage(cnd)), "134")
})

test_that("a variance condition pinned well off the truth is caught", {
  expect_warning(
    fit <- gmm_estimate(
      stacked_equations = make_mean_variance_psi(variance = 5),
      init = 2,
      overid_maxiter = settled_budget
    ),
    class = "deli_gmm_moments_rejected"
  )
  expect_lt(j_p_value(fit), 1e-8)
})

test_that("the reported understated-variance fit is caught", {
  # This system's updating loop is slow, so it is given room to settle: J is
  # only interpretable once the weight matrix has.
  psi <- make_understated_variance_psi()
  expect_warning(
    fit <- gmm_estimate(
      stacked_equations = psi,
      init = 2,
      overid_maxiter = 4000L
    ),
    class = "deli_gmm_moments_rejected"
  )
  expect_gt(max(abs(rowSums(psi(fit@theta)))), 10000)
  expect_lt(j_p_value(fit), 1e-8)
})

# ---- the threshold -----------------------------------------------------------
#
# The warning has to catch a specification failure without ever firing on an
# honest fit, and the chi-squared reference is asymptotic, so at small n and at a
# weight matrix that has not settled J is over-sized. Measured over two hundred
# correctly specified fits at each size, the smallest P-value was 1.1e-4 at
# n = 15, 2.2e-5 at n = 30 and 2.9e-6 at n = 60; at the identity weight matrix it
# was 4.8e-6. The weakest signal from a badly specified fit was 8.8e-11. The
# threshold sits in the gap between those, near its logarithmic middle, so it is
# two orders of magnitude clear of the worst honest fit measured and two orders
# clear of the weakest real failure.

test_that("a small over-identified sample does not trip the J warning", {
  # n = 15, where the chi-squared reference is at its least accurate. Seed 1011
  # of this set reports a P-value of 6.3e-4 while being correctly specified,
  # which is what rules out any threshold near a thousandth. Seed 1007 is left
  # out because its updating loop does not settle at any budget tried, so it
  # carries an over-identification warning of its own and says nothing about the
  # J threshold.
  for (seed in c(1001:1006, 1008:1012)) {
    set.seed(seed)
    y <- stats::rnorm(15, mean = 3, sd = sqrt(2))
    psi <- function(theta) rbind(y - theta[1], (y - theta[1])^2 - 2)
    expect_no_warning(gmm_estimate(
      stacked_equations = psi,
      init = 3,
      overid_maxiter = settled_budget
    ))
  }
})

test_that("an identity weight matrix does not trip the J warning", {
  # `overid_maxiter = 0` skips the updating, so J is left standing on the
  # identity, where it is an unstandardized sum of squared moments and its
  # chi-squared reference does not hold. The measured smallest P-value over a
  # hundred correctly specified fits there is 4.8e-6, which the threshold has to
  # clear as well.
  expect_no_warning(gmm_estimate(
    stacked_equations = make_mean_variance_psi(),
    init = 2,
    overid_maxiter = 0L
  ))
  expect_no_warning(gmm_estimate(
    stacked_equations = make_poisson_psi(),
    init = 1,
    overid_maxiter = 0L
  ))
})

# ---- the formula interface carries J too -------------------------------------

test_that("a formula GMM fit carries the J-statistic", {
  set.seed(13)
  n <- 400
  d <- data.frame(
    x1 = stats::rnorm(n),
    x2 = stats::rnorm(n),
    spare = stats::rnorm(n)
  )
  d$y <- 1 + 2 * d$x1 + 3 * d$x2 + stats::rnorm(n, sd = 0.5)
  # A fourth moment condition from a variable the outcome does not depend on, so
  # the model is correct and the extra condition is valid.
  ee_extra_moment <- function(theta, X, y, spare, ...) {
    residual <- as.numeric(y - X %*% theta)
    rbind(t(X * residual), spare * residual)
  }
  fit <- expect_no_warning(gmm_estimate(
    y ~ x1 + x2,
    data = d,
    .ee = ee_extra_moment,
    spare = spare,
    overid_maxiter = settled_budget
  ))
  expect_identical(j_df(fit), 1L)
  expect_true(is.finite(fit@j_statistic))
  # The measured P-value is 0.535.
  expect_gt(j_p_value(fit), 0.05)
})

# ---- moment conditions that repeat one another -------------------------------
#
# A stack whose moment conditions are linearly dependent has a singular moment
# covariance, so the two-step update cannot invert it and reaches
# `MASS::ginv()`. That happened without a word: the fit came back with the
# estimates and the standard errors of the independent conditions on their own,
# and the J-statistic said nothing either, because a repeated condition adds no
# disagreement for J to measure and drives it to zero.
#
# The reading here is about the weight matrix alone. Whether the parameters are
# identified is a separate question that the bread answers, and this stack
# identifies them perfectly well.

dependent_moment_case <- function() {
  set.seed(7)
  n <- 200
  x <- stats::rnorm(n)
  y <- 1 + 2 * x + stats::rnorm(n)
  X <- cbind(1, x)
  # The third condition repeats the first, so three conditions carry the
  # information of two and the moment covariance is rank two of three.
  function(theta) {
    r <- as.numeric(y - X %*% theta)
    rbind(r, r * x, r)
  }
}

test_that("a repeated moment condition is reported rather than pseudo-inverted", {
  skip_if_not_installed("MASS")
  psi <- dependent_moment_case()

  expect_warning(
    fit <- gmm_estimate(
      stacked_equations = psi,
      init = c(0, 0),
      overid_maxiter = settled_budget
    ),
    class = "deli_gmm_moments_dependent"
  )

  # The estimates are still the ones the objective asked for, so this is a
  # report on the system rather than a refusal to fit it.
  expect_length(coef(fit), 2L)
  expect_equal(unname(coef(fit)), c(0.9477646, 2.0693459), tolerance = 1e-6)
})

test_that("the report names the conditions that carry no new information", {
  skip_if_not_installed("MASS")
  psi <- dependent_moment_case()

  w <- expect_warning(
    gmm_estimate(
      stacked_equations = psi,
      init = c(0, 0),
      overid_maxiter = settled_budget
    ),
    class = "deli_gmm_moments_dependent"
  )

  expect_match(conditionMessage(w), "3")
  expect_match(conditionMessage(w), "rank 2 of 3")
})

test_that("the J-statistic does not surface a repeated moment condition", {
  skip_if_not_installed("MASS")
  psi <- dependent_moment_case()

  expect_warning(
    fit <- gmm_estimate(
      stacked_equations = psi,
      init = c(0, 0),
      overid_maxiter = settled_budget
    ),
    class = "deli_gmm_moments_dependent"
  )

  # A repeated condition agrees with the one it repeats by construction, so J
  # is at round-off and its P-value is 1. Nothing about the reading J takes
  # would ever have reported this stack.
  expect_lt(fit@j_statistic, 1e-8)
  expect_equal(j_p_value(fit), 1, tolerance = 1e-6)
})

test_that("an independent over-identified stack raises no dependence warning", {
  skip_if_not_installed("MASS")
  set.seed(7)
  n <- 200
  x <- stats::rnorm(n)
  y <- 1 + 2 * x + stats::rnorm(n)
  X <- cbind(1, x)
  psi <- function(theta) {
    r <- as.numeric(y - X %*% theta)
    rbind(r, r * x, r * x^2)
  }

  expect_no_warning(
    gmm_estimate(
      stacked_equations = psi,
      init = c(0, 0),
      overid_maxiter = settled_budget
    )
  )
})

# ---- shaq_free_throws --------------------------------------------------------

test_that("shaq_free_throws is a data.frame with correct structure", {
  expect_s3_class(shaq_free_throws, "data.frame")
  expect_equal(nrow(shaq_free_throws), 23)
  expect_equal(ncol(shaq_free_throws), 3)
  expect_named(shaq_free_throws, c("game", "ft_success", "ft_attempt"))
})

test_that("shaq_free_throws values are non-negative integers", {
  expect_true(all(shaq_free_throws$ft_success >= 0))
  expect_true(all(shaq_free_throws$ft_attempt >= 0))
  expect_true(all(shaq_free_throws$ft_success == as.integer(shaq_free_throws$ft_success)))
  expect_true(all(shaq_free_throws$ft_attempt == as.integer(shaq_free_throws$ft_attempt)))
})

test_that("shaq_free_throws made <= attempts for every row", {
  expect_true(all(shaq_free_throws$ft_success <= shaq_free_throws$ft_attempt))
})

test_that("shaq_free_throws first row matches Python source", {
  expect_equal(shaq_free_throws$game[1], 1)
  expect_equal(shaq_free_throws$ft_success[1], 4)
  expect_equal(shaq_free_throws$ft_attempt[1], 5)
})

# ---- inderjit ----------------------------------------------------------------

test_that("inderjit is a data.frame with correct structure", {
  expect_s3_class(inderjit, "data.frame")
  expect_equal(nrow(inderjit), 24)
  expect_equal(ncol(inderjit), 2)
  expect_named(inderjit, c("response", "dose"))
})

test_that("inderjit dose values are non-negative", {
  expect_true(all(inderjit$dose >= 0))
})

test_that("inderjit first and last rows match Python source", {
  expect_equal(inderjit$response[1], 7.58)
  expect_equal(inderjit$dose[1], 0)
  expect_equal(inderjit$response[24], 0.44)
  expect_equal(inderjit$dose[24], 30)
})

# ---- robust_regress ----------------------------------------------------------

test_that("robust_regress is a data.frame with correct structure", {

  expect_s3_class(robust_regress, "data.frame")
  expect_equal(nrow(robust_regress), 15)
  expect_equal(ncol(robust_regress), 3)
  expect_named(robust_regress, c("height", "weight", "weight_no_outlier"))
})

test_that("robust_regress values match Python source", {
  # First row

expect_equal(robust_regress$height[1], 168.519)
  expect_equal(robust_regress$weight_no_outlier[1], 67.634)

  # Row 9 has the outlier: weight = 62.043 + 3, weight_no_outlier = 62.043
  expect_equal(robust_regress$weight[9], 62.043 + 3)
  expect_equal(robust_regress$weight_no_outlier[9], 62.043)
})

# ---- breast_cancer -----------------------------------------------------------

test_that("breast_cancer is a data.frame with correct structure", {
  expect_s3_class(breast_cancer, "data.frame")
  expect_equal(nrow(breast_cancer), 45)
  expect_equal(ncol(breast_cancer), 3)
  expect_named(breast_cancer, c("delta", "times", "stain"))
})

test_that("breast_cancer binary columns contain only 0/1", {
  expect_true(all(breast_cancer$delta %in% c(0, 1)))
  expect_true(all(breast_cancer$stain %in% c(0, 1)))
})

test_that("breast_cancer first row matches Python source", {
  expect_equal(breast_cancer$delta[1], 1)
  expect_equal(breast_cancer$times[1], 23)
  expect_equal(breast_cancer$stain[1], 0)
})

# ---- sdss_quasar -------------------------------------------------------------

test_that("sdss_quasar retains the full SDSS 3rd Data Release table", {
  expect_s3_class(sdss_quasar, "data.frame")
  expect_equal(nrow(sdss_quasar), 46420)
  expect_equal(ncol(sdss_quasar), 23)
  expect_true(all(c("z", "u_mag", "g_mag", "r_mag", "i_mag") %in% names(sdss_quasar)))
})

test_that("sdss_quasar first row values are preserved", {
  expect_equal(sdss_quasar$SDSS_J[1], "000009.26+151754.5")
  expect_equal(sdss_quasar$z[1], 1.1986)
  expect_equal(sdss_quasar$u_mag[1], 19.921)
})

# ---- nsduh -------------------------------------------------------------------

test_that("nsduh retains the full survey", {
  expect_s3_class(nsduh, "data.frame")
  expect_equal(nrow(nsduh), 40293)
  expect_equal(ncol(nsduh), 14)
  expect_true(all(c("cigmon", "lgb_flag", "cig15", "age", "race") %in% names(nsduh)))
})

test_that("nsduh first row values are preserved", {
  expect_equal(nsduh$cigmon[1], 0L)
  expect_equal(nsduh$educ[1], 2L)
  expect_equal(nsduh$income[1], 4L)
})

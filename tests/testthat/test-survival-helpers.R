# ---- survival returns input unchanged ----------------------------------------
test_that("convert_survival_measures returns survival unchanged", {
  expect_equal(convert_survival_measures(0.8, measure = "survival"), 0.8)
  expect_equal(convert_survival_measures(0.5, measure = "survival"), 0.5)
  expect_equal(convert_survival_measures(1.0, measure = "survival"), 1.0)
})

# ---- risk returns 1 - survival ----------------------------------------------
test_that("convert_survival_measures returns 1 - survival for risk", {
  expect_equal(convert_survival_measures(0.8, measure = "risk"), 0.2)
  expect_equal(convert_survival_measures(0.5, measure = "risk"), 0.5)
  expect_equal(convert_survival_measures(1.0, measure = "risk"), 0.0)
})

# ---- cumulative hazard returns -log(survival) --------------------------------
test_that("convert_survival_measures returns -log(survival) for cumulative_hazard", {
  expect_equal(
    convert_survival_measures(0.8, measure = "cumulative_hazard"),
    -log(0.8)
  )
  expect_equal(
    convert_survival_measures(0.5, measure = "cumulative_hazard"),
    -log(0.5)
  )
  expect_equal(
    convert_survival_measures(1.0, measure = "cumulative_hazard"),
    0.0
  )
})

# ---- vector input works for each measure -------------------------------------
test_that("convert_survival_measures handles vector input for survival", {
  surv <- c(0.9, 0.75, 0.5, 0.25)
  expect_equal(
    convert_survival_measures(surv, measure = "survival"),
    surv
  )
})

test_that("convert_survival_measures handles vector input for risk", {
  surv <- c(0.9, 0.75, 0.5, 0.25)
  expect_equal(
    convert_survival_measures(surv, measure = "risk"),
    1 - surv
  )
})

test_that("convert_survival_measures handles vector input for cumulative_hazard", {
  surv <- c(0.9, 0.75, 0.5, 0.25)
  expect_equal(
    convert_survival_measures(surv, measure = "cumulative_hazard"),
    -log(surv)
  )
})

# ---- invalid measure raises an error -----------------------------------------
test_that("convert_survival_measures errors on invalid measure", {
  expect_error(
    convert_survival_measures(0.8, measure = "unknown"),
    "measure"
  )
})

# ---- edge cases --------------------------------------------------------------
test_that("convert_survival_measures handles survival = 1 edge case", {
  expect_equal(convert_survival_measures(1.0, measure = "risk"), 0.0)
  expect_equal(convert_survival_measures(1.0, measure = "cumulative_hazard"), 0.0)
})

test_that("convert_survival_measures handles survival near 0", {
  near_zero <- 1e-10
  expect_equal(
    convert_survival_measures(near_zero, measure = "risk"),
    1 - near_zero
  )
  expect_equal(
    convert_survival_measures(near_zero, measure = "cumulative_hazard"),
    -log(near_zero)
  )
  # cumulative hazard should be a large positive number when survival is near 0
  expect_true(
    convert_survival_measures(near_zero, measure = "cumulative_hazard") > 0
  )
})

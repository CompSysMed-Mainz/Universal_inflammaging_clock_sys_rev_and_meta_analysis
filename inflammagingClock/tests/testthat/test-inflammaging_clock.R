make_test_clock <- function() {
  list(
    features = c("IL6", "TNF", "CXCL10"),
    training_medians = c(IL6 = 1, TNF = 5, CXCL10 = 3),
    max_absent_feature_fraction = 0.4,
    model_type = "elastic_net",
    intercept = 20,
    coefficients = c(CXCL10 = -1, TNF = 0.5, IL6 = 2),
    age_calibration = c(2, 0.8)
  )
}

test_that("duplicate gene assays are averaged and coefficients matched", {
  samples <- data.frame(
    id = c("A", "B"), age = c(70, 60),
    IL6 = c(1, 4), ` il6 ` = c(3, 6),
    TNF = c(NA_real_, 4), CXCL10 = c(2, 7),
    check.names = FALSE
  )
  scores <- inflammaging_clock(samples, make_test_clock(), "id", "age")

  expect_equal(scores$predicted_age, c(24.5, 25))
  expect_equal(scores$imputed_features, c(1, 0))
  expect_equal(scores$absent_feature_fraction, c(0, 0))
  expect_equal(scores$age_acceleration, c(24.5 - 58, 25 - 50))
  expect_equal(scores$sample_id, c("A", "B"))
})

test_that("absent features use training medians", {
  samples <- data.frame(IL6 = c(1, 2), TNF = c(2, 3))
  scores <- inflammaging_clock(samples, make_test_clock())

  expect_equal(scores$absent_feature_fraction, rep(1 / 3, 2))
  expect_equal(scores$imputed_features, c(1, 1))
  expect_equal(scores$predicted_age, c(20, 22.5))
})

test_that("too many absent features are rejected", {
  expect_error(
    inflammaging_clock(data.frame(IL6 = 1), make_test_clock()),
    "absent"
  )
})

test_that("assay and age columns must be numeric", {
  expect_error(
    inflammaging_clock(data.frame(IL6 = 1, TNF = "bad", CXCL10 = 2),
                       make_test_clock()),
    "numeric"
  )
  expect_error(
    inflammaging_clock(data.frame(age = factor(70), IL6 = 1,
                                   TNF = 2, CXCL10 = 3),
                       make_test_clock(), age_column = "age"),
    "age column must be numeric"
  )
})

test_that("missing coefficients are rejected", {
  clock <- make_test_clock()
  clock$coefficients <- c(IL6 = 2, TNF = 0.5)
  expect_error(
    inflammaging_clock(data.frame(IL6 = 1, TNF = 2, CXCL10 = 3), clock),
    "named coefficients"
  )
})

test_that("the bundled RDS can be loaded", {
  model_path <- system.file(
    "extdata", "inflammaging_clock.rds", package = "inflammagingClock"
  )
  skip_if(!nzchar(model_path), "The validated RDS has not been added yet")

  clock <- load_inflammaging_clock()
  expect_true(is.list(clock))
  expect_true(length(clock$features) > 0L)
})

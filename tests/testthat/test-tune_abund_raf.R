require(dplyr)
data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

raf_grid <- expand.grid(
  mtry = seq(from = 2, to = 3, by = 1),
  ntree = seq(from = 500, to = 1000, by = 300)
)


test_that("tune_abund_raf", {
  set.seed(123)
  tuned_raf <- tune_abund_raf(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = raf_grid,
    n_cores = 1
  )

  dim(tuned_raf$optimal_combination) %>% expect_equal(c(1, 16))
  expect_true(tuned_raf$performance$corr_spear_mean < 0.5)
})

test_that("test errors", {
  expect_error(tune_abund_raf(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear", "mae"),
    grid = raf_grid,
    n_cores = 1
  ))

  expect_error(tune_abund_raf(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = raf_grid <- expand.grid(
      mtryE = seq(from = 2, to = 3, by = 1),
      ntreeE = c(100, 500)
    ),
    n_cores = 1
  ))
})

test_that("incomplete gird", {
  tuned_raf <- tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtry = seq(from = 2, to = 3, by = 1)
    ),
    n_cores = 1,
    verbose = FALSE
  )

  length(tuned_raf) %>% expect_equal(7)
})

test_that("message", {
  expect_message(tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtry = seq(from = 2, to = 3, by = 1)
    ),
    n_cores = 1,
    verbose = FALSE
  ))
})

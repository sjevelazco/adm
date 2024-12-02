require(dplyr)
data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
# Create a grid
xgb_grid <- expand.grid(
  nrounds = c(100),
  max_depth = c(8),
  eta = c(0.2, 0.5),
  gamma = c(1, 10),
  colsample_bytree = c(0.5),
  min_child_weight = c(1),
  subsample = c(0.5, 1)
)

test_that("tune_abund_xgb and fit_abund_xgb", {
  set.seed(1)
  tuned_ <- tune_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = xgb_grid,
    objective = "reg:squarederror",
    n_cores = 3
  )
  expect_equal(names(tuned_), c(
    "model", "predictors", "performance", "performance_part",
    "predicted_part", "optimal_combination", "all_combinations"
  ))
  expect_equal(class(tuned_$model)[1], "xgb.Booster")
  expect_equal(round(tuned_$performance$corr_spear_mean, 2), 0.18)
  expect_equal(tuned_$optimal_combination$nrounds, 100)
})

test_that("test errors", {
  expect_error(tune_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear", "mae"),
    grid = xgb_grid,
    objective = "reg:squarederror",
    n_cores = 3
  ))
  
  expect_error(tune_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtryE = seq(from = 2, to = 3, by = 1),
      ntreeE = c(100, 500)
    ),
    objective = "reg:squarederror",
    n_cores = 3
  ))
})

test_that("incomplete grid", {
  tuned_ <- tune_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      nrounds = c(100),
      max_depth = c(8),
      eta = c(0.2, 0.5),
      gamma = c(1, 10),
      colsample_bytree = c(0.5),
      min_child_weight = c(1)
      # subsample = c(0.5, 1)
    ),
    objective = "reg:squarederror",
    n_cores = 3
  )
  
  expect_true("subsample" %in% names(tuned_$optimal_combination))
})


require(dplyr)
data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>% 
  dplyr::select(-.part2, -.part3)
svm_grid <- expand.grid(
  sigma = "automatic",
  C = c(0.5, 2),
  kernel = c("rbfdot","laplacedot")
)

test_that("tune_abund_svm and fit_aund_svm", {
  
  set.seed(1)
  tuned_svm <- tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear","mae"),
    grid = svm_grid,
    n_cores = 1
  )
  expect_equal(names(tuned_svm), c("model", "predictors", "performance", "performance_part", 
                                   "predicted_part", "optimal_combination", "all_combinations"))
  expect_equal(class(tuned_svm$model)[1], "ksvm")
  expect_equal(round(tuned_svm$performance$corr_spear_mean, 2), 0.58)
  expect_equal(tuned_svm$optimal_combination$C, 0.5)
})

test_that("test errors", {
  expect_error(tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear","mae"),
    grid = svm_grid,
    n_cores = 1
  ))
  
  expect_message(tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear","mae"),
    grid = svm_grid[, 1:2],
    n_cores = 1,
    verbose = TRUE
  ))
 
})

test_that("message", {
  expect_message(tune_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear","mae"),
    grid = svm_grid,
    n_cores = 1,
    verbose = TRUE
  ))
})


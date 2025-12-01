require(dplyr)

# install torch

data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)
# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Architecture

test_that("tune_abund_dnn and fit_abund_dnn", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  one_arch <- generate_dnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    number_of_hidden_layers = 3,
    hidden_layers_size = c(8, 16, 8),
    batch_norm = TRUE
  )

  # Create a grid
  dnn_grid <- expand.grid(
    learning_rate = c(0.01),
    n_epochs = c(50),
    batch_size = c(32),
    validation_patience = c(2, 4),
    fitting_patience = c(2, 4)
  )

  set.seed(1)
  tuned_ <- tune_abund_dnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = dnn_grid,
    architectures = one_arch,
    n_cores = 1, verbose = FALSE
  )
  expect_equal(names(tuned_), c(
    "model", "predictors", "performance", "performance_part",
    "predicted_part", "optimal_combination", "all_combinations", "selected_arch",
    "monitor"
  ))
  expect_equal(class(tuned_$model)[1], "luz_module_fitted")
})

test_that("test errors", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  one_arch <- generate_dnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    number_of_hidden_layers = 3,
    hidden_layers_size = c(8, 16, 8),
    batch_norm = TRUE
  )

  # Create a grid
  dnn_grid <- expand.grid(
    learning_rate = c(0.01),
    n_epochs = c(50),
    batch_size = c(32),
    validation_patience = c(2, 4),
    fitting_patience = c(2, 4)
  )

  expect_error(tune_abund_dnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear", "mae"),
    grid = dnn_grid,
    architectures = one_arch,
    n_cores = 1
  ))

  expect_error(tune_abund_dnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtryE = seq(from = 2, to = 3, by = 1),
      ntreeE = c(100, 500)
    ),
    architectures = one_arch,
    n_cores = 1
  ))
})

# test_that("incomplete grid", {
#   if(!torch::torch_is_installed()){
#     skip()
#   }
#   one_arch <- generate_dnn_architecture(
#     number_of_features = 3,
#     number_of_outputs = 1,
#     number_of_hidden_layers = 3,
#     hidden_layers_size = c(8, 16, 8),
#     batch_norm = TRUE
#   )
#
#
#   if(!torch::torch_is_installed()){
#     skip()
#   }
#   set.seed(1)
#   tuned_ <- tune_abund_dnn(
#     data = some_sp,
#     response = "ind_ha",
#     predictors = c("bio12", "elevation", "sand"),
#     partition = ".part",
#     predict_part = TRUE,
#     metrics = c("corr_pear", "mae"),
#     grid = expand.grid(
#       learning_rate = c(0.01),
#       # n_epochs = c(50),
#       batch_size = c(32),
#       validation_patience = c(2,4),
#       fitting_patience = c(2,4)
#     )
#     ,
#     architectures = one_arch,
#     n_cores = 1
#   )
#
#   expect_true("n_epochs" %in% names(tuned_$optimal_combination))
#
# })

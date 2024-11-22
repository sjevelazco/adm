require(dplyr)

#install torch
torch::install_torch()

data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

one_arch <- generate_cnn_architecture(
  number_of_features = 3,
  number_of_outputs = 1,
  sample_size = c(11, 11),
  number_of_conv_layers = 2,
  conv_layers_size = c(14, 28),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  number_of_fc_layers = 1,
  fc_layers_size = c(28),
  pooling = NULL,
  batch_norm = TRUE,
  dropout = 0,
  verbose = T
)

# Create a grid
# Obs.: the grid is tested with every architecture, thus it can get very large.
cnn_grid <- expand.grid(
  learning_rate = c(0.01),
  n_epochs = c(50),
  batch_size = c(32),
  validation_patience = c(4),
  fitting_patience = c(4)
)


test_that("tune_abund_svm and fit_abund_svm", {
  set.seed(1)
  tuned_ <- tune_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = cnn_grid,
    rasters = system.file("external/envar.tif", package = "adm"),
    x = "x",
    y = "y",
    sample_size = c(11,11),
    architectures = one_arch,
    n_cores = 3,
    verbose = FALSE
  )
  expect_equal(names(tuned_), c(
    "model", "predictors", "performance", "performance_part",
    "predicted_part", "optimal_combination", "all_combinations", "selected_arch"
  ))
  expect_equal(class(tuned_$model)[1], "luz_module_fitted")
})

test_that("test errors", {
  expect_error( tune_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear", "mae"),
    grid = cnn_grid,
    rasters = system.file("external/envar.tif", package = "adm"),
    x = "x",
    y = "y",
    sample_size = c(11,11),
    architectures = one_arch,
    n_cores = 3,
    verbose = FALSE
  ))
  
  expect_error(tune_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtryE = seq(from = 2, to = 3, by = 1),
      ntreeE = c(100, 500)
    ),
    rasters = system.file("external/envar.tif", package = "adm"),
    x = "x",
    y = "y",
    sample_size = c(11,11),
    architectures = one_arch,
    n_cores = 3,
    verbose = FALSE
  ))
})

test_that("incomplete grid", {
  tuned_ <- tune_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      learning_rate = c(0.01),
      n_epochs = c(50),
      batch_size = c(32),
      validation_patience = c(4)
      # fitting_patience = c(4)
    ),
    rasters = system.file("external/envar.tif", package = "adm"),
    x = "x",
    y = "y",
    sample_size = c(11,11),
    architectures = one_arch,
    n_cores = 3,
    verbose = FALSE
  )
  
  expect_true("fitting_patience" %in% names(tuned_$optimal_combination))
  
})

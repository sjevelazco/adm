require(dplyr)

# install torch

data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
# Small subset to keep tests fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 200)), ]




test_that("tune_abund_cnn and fit_abund_cnn", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  one_arch <- generate_cnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    sample_size = c(5, 5),
    number_of_conv_layers = 2,
    conv_layers_size = c(4, 8),
    conv_layers_kernel = 3,
    conv_layers_stride = 1,
    conv_layers_padding = 0,
    number_of_fc_layers = 1,
    fc_layers_size = c(8),
    pooling = NULL,
    batch_norm = FALSE,
    dropout = 0,
    verbose = FALSE
  )

  # Create a grid
  # Obs.: the grid is tested with every architecture, thus it can get very large.
  cnn_grid <- expand.grid(
    learning_rate = c(0.0001),
    n_epochs = c(2),
    batch_size = c(32),
    validation_patience = c(2),
    fitting_patience = c(2)
  )

  set.seed(1)
  torch::torch_manual_seed(1)
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
    sample_size = c(5, 5),
    architectures = one_arch,
    n_cores = 1,
    verbose = FALSE
  )
  expect_equal(names(tuned_), c(
    "model", "predictors", "performance", "performance_part",
    "predicted_part", "metadata", "optimal_combination", "all_combinations", "selected_arch"
  ))
  expect_equal(class(tuned_$model)[1], "luz_module_fitted")
})

test_that("test errors", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  one_arch <- generate_cnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    sample_size = c(5, 5),
    number_of_conv_layers = 2,
    conv_layers_size = c(4, 8),
    conv_layers_kernel = 3,
    conv_layers_stride = 1,
    conv_layers_padding = 0,
    number_of_fc_layers = 1,
    fc_layers_size = c(8),
    pooling = NULL,
    batch_norm = FALSE,
    dropout = 0,
    verbose = FALSE
  )

  # Create a grid
  # Obs.: the grid is tested with every architecture, thus it can get very large.
  cnn_grid <- expand.grid(
    learning_rate = c(0.0001),
    n_epochs = c(2),
    batch_size = c(32),
    validation_patience = c(2),
    fitting_patience = c(2)
  )

  expect_error(tune_abund_cnn(
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
    sample_size = c(5, 5),
    architectures = one_arch,
    n_cores = 1,
    verbose = FALSE
  ))

  expect_error(
    tune_abund_cnn(
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
      sample_size = c(5, 5),
      architectures = one_arch,
      n_cores = 1,
      verbose = FALSE
    )
  )
})

test_that("incomplete grid", {
  skip_on_cran()
  if (!torch::torch_is_installed()) {
    skip()
  }

  one_arch <- generate_cnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    sample_size = c(5, 5),
    number_of_conv_layers = 2,
    conv_layers_size = c(4, 8),
    conv_layers_kernel = 3,
    conv_layers_stride = 1,
    conv_layers_padding = 0,
    number_of_fc_layers = 1,
    fc_layers_size = c(8),
    pooling = NULL,
    batch_norm = FALSE,
    dropout = 0,
    verbose = FALSE
  )

  # Create a grid
  # Obs.: the grid is tested with every architecture, thus it can get very large.
  cnn_grid <- expand.grid(
    learning_rate = c(0.0001),
    n_epochs = c(2),
    batch_size = c(32),
    validation_patience = c(2),
    fitting_patience = c(2)
  )


  tuned_ <- tune_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      learning_rate = c(0.0001),
      n_epochs = c(2),
      batch_size = c(32),
      validation_patience = c(2)
      # fitting_patience = c(2)
    ),
    rasters = system.file("external/envar.tif", package = "adm"),
    x = "x",
    y = "y",
    sample_size = c(5, 5),
    architectures = one_arch,
    n_cores = 1,
    verbose = FALSE
  )

  expect_true("fitting_patience" %in% names(tuned_$optimal_combination))
})

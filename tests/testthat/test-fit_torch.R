# Quick DNN/CNN fits (tiny data and epochs); skipped when libtorch is missing
set.seed(1)
torch_sp <- test_sp[sort(sample(nrow(test_sp), 120)), ]

test_that("fit_abund_dnn default architecture, partition and hold-out", {
  skip_if_not(torch::torch_is_installed())
  m <- suppressMessages(fit_abund_dnn(torch_sp, "ind_ha", test_pred,
    partition = ".part1", predict_part = TRUE, hold_out_set = torch_sp[1:20, ],
    n_epochs = 2, validation_patience = 1, fitting_patience = 1, verbose = FALSE
  ))
  expect_s3_class(m$model, "luz_module_fitted")
  expect_true(all(c("performance", "performance_part", "performance_part_hold_out") %in% names(m)))
  expect_false(is.null(m$predicted_part))
})

test_that("fit_abund_dnn without partition, warning for categorical predictors, monitor", {
  skip_if_not(torch::torch_is_installed())
  expect_warning(
    m <- suppressMessages(fit_abund_dnn(torch_sp, "ind_ha", test_pred,
      predictors_f = "eco", partition = NULL,
      n_epochs = 2, fitting_patience = 1, verbose = FALSE
    )),
    "Categorical"
  )
  expect_true(all(c("model", "predictors", "metadata") %in% names(m)))

  m <- suppressMessages(fit_abund_dnn(torch_sp, "ind_ha", test_pred,
    partition = ".part1", n_epochs = 2, validation_patience = 1,
    fitting_patience = 1, learning_monitor = TRUE, verbose = FALSE
  ))
  expect_s3_class(m$model, "luz_module_fitted")
})

test_that("fit_abund_cnn without partition; hold-out is rejected", {
  skip_on_cran()
  skip_if_not(torch::torch_is_installed())
  arch <- generate_cnn_architecture(
    number_of_features = 3, number_of_outputs = 1, sample_size = c(5, 5),
    number_of_conv_layers = 1, conv_layers_size = 4, conv_layers_kernel = 3,
    number_of_fc_layers = 1, fc_layers_size = 4, batch_norm = FALSE
  )
  m0 <- suppressMessages(fit_abund_cnn(torch_sp, "ind_ha", test_pred,
    partition = NULL, x = "x", y = "y", rasters = test_envar,
    sample_size = c(5, 5), n_epochs = 2, fitting_patience = 1,
    custom_architecture = arch, verbose = FALSE
  ))
  expect_true(all(c("model", "predictors", "metadata") %in% names(m0)))

  m <- suppressWarnings(suppressMessages(fit_abund_cnn(torch_sp, "ind_ha", test_pred,
    partition = ".part1", x = "x", y = "y", rasters = test_envar,
    sample_size = c(5, 5), n_epochs = 2, validation_patience = 1,
    fitting_patience = 1, custom_architecture = arch, verbose = FALSE
  )))
  expect_s3_class(m$model, "luz_module_fitted")

  expect_error(
    fit_abund_cnn(torch_sp, "ind_ha", test_pred,
      partition = ".part1", x = "x", y = "y", rasters = test_envar,
      sample_size = c(5, 5), hold_out_set = torch_sp[1:20, ], verbose = FALSE
    ),
    "hold_out_set"
  )
})

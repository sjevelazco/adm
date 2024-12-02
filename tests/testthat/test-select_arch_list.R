# install torch
if (!torch::torch_is_installed()) {
  skip()
}

test_that("select_arch_list for dnn", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  archs <- generate_arch_list(
    type = "dnn",
    number_of_features = 4,
    number_of_outputs = 1,
    n_layers = c(2, 3, 4),
    n_neurons = c(8, 16, 64),
    batch_norm = TRUE,
    dropout = 0.2
  ) %>% select_arch_list(
    type = c("dnn"),
    method = "percentile",
    n_samples = 1,
    min_max = TRUE # Keep the network with the minimum and maximum number of parameters
  )

  expect_equal(names(archs), c("arch_list", "arch_dict", "changes"))
  expect_equal(length(archs$arch_dict), 3)
  expect_equal(length(archs$arch_list), 28)
})

test_that("select_arch_list for dnn, min_max FALSE", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  archs <- generate_arch_list(
    type = "dnn",
    number_of_features = 4,
    number_of_outputs = 1,
    n_layers = c(2, 3, 4),
    n_neurons = c(8, 16, 64),
    batch_norm = TRUE,
    dropout = 0.2
  ) %>% select_arch_list(
    type = c("dnn"),
    method = "percentile",
    n_samples = 1,
    min_max = FALSE # Keep the network with the minimum and maximum number of parameters
  )

  expect_equal(names(archs), c("arch_list", "arch_dict", "changes"))
  expect_equal(length(archs$arch_dict), 3)
  expect_equal(length(archs$arch_list), 26)
})

test_that("select_arch_list for cnn", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  archs <- generate_arch_list(
    type = "cnn",
    number_of_features = 4,
    number_of_outputs = 1,
    n_layers = c(2, 4), # now convolutional layers
    n_neurons = c(8, 64),
    sample_size = c(11, 11),
    number_of_fc_layers = c(2, 4), # fully connected layers
    fc_layers_size = c(16, 8),
    conv_layers_kernel = 3,
    conv_layers_stride = 1,
    conv_layers_padding = 0,
    pooling = 1,
    batch_norm = TRUE,
    dropout = 0.2
  ) %>% select_arch_list(
    type = c("cnn"),
    method = "percentile",
    n_samples = 1,
    min_max = TRUE # Keep the network with the minimum and maximum number of parameters
  )

  expect_equal(names(archs), c("arch_list", "arch_dict", "changes"))
  expect_equal(length(archs$arch_dict), 4)
  expect_equal(length(archs$arch_list), 38)
})

test_that("select_arch_list for cnn dropout", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  archs <- generate_arch_list(
    type = "cnn",
    number_of_features = 4,
    number_of_outputs = 1,
    n_layers = c(2, 4), # now convolutional layers
    n_neurons = c(8),
    sample_size = c(11, 11),
    number_of_fc_layers = c(2, 4), # fully connected layers
    fc_layers_size = c(8),
    conv_layers_kernel = 3,
    conv_layers_stride = 1,
    conv_layers_padding = 0,
    pooling = 1,
    batch_norm = TRUE,
    dropout = 0.2
  ) %>% select_arch_list(
    type = c("cnn"),
    method = "percentile",
    n_samples = 1,
    min_max = TRUE # Keep the network with the minimum and maximum number of parameters
  )
  expect_equal(names(archs), c("arch_list", "arch_dict", "changes"))
  expect_equal(length(archs$arch_dict), 4)
  expect_equal(length(archs$arch_list), 4)
})

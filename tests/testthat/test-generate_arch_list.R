#install torch
torch::install_torch()

test_that("architecture for cnn", {
  dnn_archs <- generate_arch_list(
    type = "dnn",
    number_of_features = 4,
    number_of_outputs = 1,
    n_layers = c(2, 3, 4),
    n_neurons = c(8, 16, 32, 64),
    batch_norm = TRUE,
    dropout = 0.2
  )

  expect_equal(names(dnn_archs), c("arch_list", "arch_dict"))
  expect_equal(as.vector(sapply(dnn_archs$arch_dict, nrow)), c(2, 3, 4))
  expect_equal(length(dnn_archs$arch_list), 336)
})


test_that("architecture for cnn", {
  cnn_archs <- generate_arch_list(
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
  )

  expect_equal(names(cnn_archs), c("arch_list", "arch_dict"))
  expect_equal(as.vector(sapply(cnn_archs$arch_dict, nrow)), c(4, 6, 6, 8))
  expect_equal(length(cnn_archs$arch_list), 400)
})

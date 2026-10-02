test_that("generate_dnn_architecture", {
  skip_if_not(torch::torch_is_installed())
  a <- generate_dnn_architecture()
  expect_true(all(c("arch", "arch_dict", "net") %in% names(a)))
  expect_true(grepl("nn_linear", a$arch))
})

test_that("generate_cnn_architecture", {
  skip_if_not(torch::torch_is_installed())
  a <- generate_cnn_architecture(sample_size = c(11, 11))
  expect_type(a, "list")
  expect_error(
    generate_cnn_architecture(sample_size = c(2, 2)),
    "too small"
  )
})

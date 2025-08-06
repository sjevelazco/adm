## %######################################################%##
#                                                          #
####      Test data_abund_bpdp and data_abund_pdp       ####
####             functions and plot verions             ####
#                                                          #
## %######################################################%##
require(dplyr)
require(terra)

# Load data
envar <- system.file("external/envar.tif", package = "adm") %>%
  rast()
data("sppabund")

some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# TODO activate this test for cnn
# test_that("data_abund_bpdp for cnn", {
#   if (!torch::torch_is_installed()) {
#     skip()
#   }
#
#   cnn_arch <- generate_cnn_architecture(
#       number_of_features = 3,
#       number_of_outputs = 1,
#       sample_size = c(11, 11),
#       number_of_conv_layers = 2,
#       conv_layers_size = c(14, 28),
#       conv_layers_kernel = 3,
#       conv_layers_stride = 1,
#       conv_layers_padding = 0,
#       number_of_fc_layers = 1,
#       fc_layers_size = c(28),
#       pooling = NULL,
#       batch_norm = TRUE,
#       dropout = 0,
#       verbose = T
#     )
#
#  set.seed(1)
#  suppressMessages(
#    m <- fit_abund_cnn(
#      data = some_sp,
#      response = "ind_ha",
#      predictors = c("bio12", "elevation", "sand"),
#      predictors_f = NULL,
#      partition = ".part",
#      x = "x",
#      y = "y",
#      rasters = envar,
#      sample_size = c(11, 11),
#      learning_rate = 0.01,
#      n_epochs = 100,
#      batch_size = 32,
#      validation_patience = 2,
#      fitting_patience = 5,
#      custom_architecture = cnn_arch,
#      verbose = TRUE,
#      predict_part = TRUE
#    )
#  )
#
# # BPDP
#   bpdp_data <- data_abund_bpdp(
#     model = m,
#     predictors = c("bio12", "sand"),
#     resolution = 25,
#     training_data = some_sp,
#     response_name = "Abundance",
#     projection_data = envar,
#     training_boundaries = "convexh"
#   )
#
#   expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
#   expect_equal(length(bpdp_data$pdpdata), 3)
#   expect_equal(ncol(bpdp_data$training_boundaries), 2)
#   expect_equal(nrow(bpdp_data$pdpdata), 25*25)
#
# # PDP
# pdp_data <- data_abund_pdp(
#   model = m,
#   predictors = "bio12",
#   resolution = 25,
#   resid = TRUE,
#   training_data = some_sp,
#   response_name = "Abundance",
#   projection_data = envar
# )
# expect_equal(names(pdp_data), c("pdpdata", "resid"))
# expect_equal(length(pdp_data$pdpdata), 3)
# expect_equal(ncol(pdp_data$resid), 2)
# expect_equal(nrow(pdp_data$pdpdata), 27)
# })

test_that("data_abund_bpdp for dnn", {
  if (!torch::torch_is_installed()) {
    skip()
  }

  # Generate a architecture
  dnn_arch <- generate_dnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    number_of_hidden_layers = 3,
    hidden_layers_size = c(8, 16, 8),
    batch_norm = TRUE
  )
  set.seed(1)
  suppressMessages(
    m <- fit_abund_dnn(
      data = some_sp,
      response = "ind_ha",
      predictors = c("bio12", "elevation", "sand"),
      predictors_f = NULL,
      partition = ".part",
      learning_rate = 0.01,
      n_epochs = 10,
      batch_size = 32,
      validation_patience = 2,
      fitting_patience = 5,
      custom_architecture = dnn_arch,
      verbose = TRUE,
      predict_part = TRUE
    )
  )
  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for gam", {
  set.seed(1)
  suppressMessages(
    m <- fit_abund_gam(
      data = some_sp,
      response = "ind_ha",
      predictors = c("elevation", "sand", "bio3", "bio12"),
      sigma_formula = ~ elevation + bio3 + bio12,
      predictors_f = NULL,
      partition = ".part",
      distribution = gamlss.dist::ZAGA()
    )
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for glm", {
  m <- fit_abund_glm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    distribution = "ZAIG",
    poly = 3,
    inter_order = 0,
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for gbm", {
  m <- fit_abund_gbm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    distribution = "gaussian",
    n.trees = 100,
    interaction.depth = 5,
    n.minobsinnode = 5,
    shrinkage = 0.1,
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for net", {
  m <- fit_abund_net(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    size = 32,
    decay = 0.1,
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for raf", {
  m <- fit_abund_raf(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    mtry = 3,
    ntree = 500,
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for svm", {
  m <- fit_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    kernel = "rbfdot",
    sigma = "automatic",
    C = 1,
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})

test_that("data_abund_bpdp for xgb", {
  m <- fit_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = NULL,
    partition = ".part",
    nrounds = 200,
    max_depth = 5,
    eta = 0.1,
    gamma = 1,
    colsample_bytree = 0.7,
    min_child_weight = 2,
    subsample = 0.3,
    objective = "reg:squarederror",
    predict_part = TRUE
  )

  # BPDP
  bpdp_data <- data_abund_bpdp(
    model = m,
    predictors = c("bio12", "sand"),
    resolution = 25,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar,
    training_boundaries = "convexh"
  )

  expect_equal(names(bpdp_data), c("pdpdata", "training_boundaries"))
  expect_equal(length(bpdp_data$pdpdata), 3)
  expect_equal(ncol(bpdp_data$training_boundaries), 2)
  expect_equal(nrow(bpdp_data$pdpdata), 25 * 25)

  # PDP
  pdp_data <- data_abund_pdp(
    model = m,
    predictors = "bio12",
    resolution = 25,
    resid = TRUE,
    training_data = some_sp,
    response_name = "Abundance",
    projection_data = envar
  )
  expect_equal(names(pdp_data), c("pdpdata", "resid"))
  expect_equal(length(pdp_data$pdpdata), 3)
  expect_equal(ncol(pdp_data$resid), 2)
  expect_equal(nrow(pdp_data$pdpdata), 27)
})


test_that("expect error", {
  m <- fit_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    kernel = "rbfdot",
    sigma = "automatic",
    C = 1,
    predict_part = TRUE
  )

  expect_error(
    data_abund_bpdp(
      model = m,
      # predictors = c("bio12", "sand"),
      resolution = 25,
      training_data = some_sp,
      response_name = "Abundance",
      projection_data = envar,
      training_boundaries = "convexh"
    )
  )

  expect_error(
    data_abund_bpdp(
      model = m,
      predictors = c("bio12", "elevation", ),
      resolution = 25,
      # training_data = some_sp,
      response_name = "Abundance",
      projection_data = envar,
      training_boundaries = "convexh"
    )
  )

  expect_error(
    data_abund_pdp(
      model = m,
      # predictors = c("bio12"),
      resolution = 25,
      training_data = some_sp,
      response_name = "Abundance",
      projection_data = envar
    )
  )

  expect_error(
    data_abund_pdp(
      model = m,
      predictors = "elevation",
      resolution = 25,
      # training_data = some_sp,
      response_name = "Abundance",
      projection_data = envar,
      resid = TRUE
    )
  )
})

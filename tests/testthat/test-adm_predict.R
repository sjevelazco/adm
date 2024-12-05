##%######################################################%##
#                                                          #
####                    Predict test                    ####
#                                                          #
##%######################################################%##
require(dplyr)
require(terra)
envar <- system.file("external/envar.tif", package = "adm") %>%
  rast()

data("sppabund")
# Extract data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

envar <- terra::crop(envar, terra::ext(as.matrix(some_sp[c("x", "y")])))
plot(envar)

#### predict XGB ####
test_that("predic XGB", {
  set.seed(123)
  m <- fit_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    # predictors_f = c("eco"),
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
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "xgb")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict SVM ####
test_that("predic SVM", {
  set.seed(123)
  m <- fit_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    kernel = "rbfdot",
    sigma = "automatic",
    C = 1,
    predict_part = TRUE
  )
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "svm")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict NET ####
test_that("predic NET", {
  set.seed(123)
  m <- fit_abund_net(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    size = 32,
    decay = 0.01,
    predict_part = TRUE
  )
  
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "net")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict GBM ####
test_that("predic GBM", {
  set.seed(123)
  m <- fit_abund_gbm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    distribution = "gaussian",
    n.trees = 100,
    interaction.depth = 5,
    n.minobsinnode = 5,
    shrinkage = 0.1,
    predict_part = TRUE
  )
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "gbm")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict GLM ####
test_that("predic GLM", {
  set.seed(123)
  m <- fit_abund_glm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = NULL,
    partition = ".part",
    distribution = "NO",
    poly = 3,
    inter_order = 2,
    predict_part = TRUE
  )
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "glm")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict GAM ####
test_that("predic GAM", {
  set.seed(123)
  suppressWarnings(m <- fit_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    distribution = "NO",
    inter = "automatic",
    predict_part = TRUE
  ))
  
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "gam")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict RAF ####
test_that("predic RAF", {
  set.seed(123)
  m <- fit_abund_raf(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
    predictors_f = c("eco"),
    partition = ".part",
    mtry = 3,
    ntree = 500,
    predict_part = TRUE
  )
  
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "raf")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict CNN ####
test_that("predic CNN", {
  if (!torch::torch_is_installed()) {
    skip()
  }
  set.seed(123)
  # Generate an architecture
  cnn_arch <- generate_cnn_architecture(
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
  
  if (!torch::torch_is_installed()) {
    skip()
  }
  # Fit a CNN model
  m <- fit_abund_cnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = NULL,
    partition = ".part",
    x = "x",
    y = "y",
    rasters = envar,
    sample_size = c(11, 11),
    learning_rate = 0.01,
    n_epochs = 100,
    batch_size = 32,
    validation_patience = 2,
    fitting_patience = 5,
    custom_architecture = cnn_arch,
    verbose = TRUE,
    predict_part = TRUE
  )
  
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "cnn")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

#### predict DNN ####
test_that("predic DNN", {
  if (!torch::torch_is_installed()) {
    skip()
  }
  set.seed(123)
  # Generate a architecture
  dnn_arch <- generate_dnn_architecture(
    number_of_features = 3,
    number_of_outputs = 1,
    number_of_hidden_layers = 3,
    hidden_layers_size = c(8, 16, 8),
    batch_norm = TRUE
  )
  
  if (!torch::torch_is_installed()) {
    skip()
  }
  # Fit a DNN model
  m <- fit_abund_dnn(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12","elevation","sand"),
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
  
  prd <- adm_predict(
    m,
    envar,
    training_data = some_sp,
    transform_negative = TRUE)
  
  expect_equal(names(prd) , "dnn")
  expect_equal(class(prd) , "list")
  expect_equal(class(prd[[1]])[[1]] , "SpatRaster")
})

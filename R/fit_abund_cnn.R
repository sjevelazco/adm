###
#' Fit and train Convolutional Neural Network Model
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param longitude character. The name of the column containing longitude information for each observation.
#' @param latitude character. The name of the column containing latitude information for each observation.
#' @param rasters a terra SpatRaster object. A raster containing the predictor variables to be cropped around each observation.
#' @param crop_size numeric. An integer defining the range of pixels around the observation from the raster object passed to rasters parameter. Default = 5
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param learning_rate numeric. The size of the step taken during the optimization process. Default = 0.01
#' @param n_epochs numeric. How many times the learning algorithm will work through the training set. Default = 10
#' @param batch_size numeric. A batch is a subset of the training set used in a single iteration of the training process. The size of each batch is referred to as the batch size. Default = 32
#'
#' @return
#' @export
#'
#' @examples
fit_abund_cnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           longitude,
           latitude,
           rasters,
           crop_size = 5,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           learning_rate = 0.01,
           n_epochs = 10,
           batch_size = 32) {
    # Variables
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()

    create_dataset <- torch::dataset(
      "dataset",
      initialize = function(data_list) {
        self$predictors <- data_list$predictors
        self$response_variable <- data_list$response
      },
      .getitem = function(index) {
        response <- torch_tensor(self$response_variable[[index]])
        x <- transform_to_tensor(self$predictors[[index]])
        list(x = x, y = response)
      },
      .length = function() {
        length(self$response_variable)
      }
    )

    ##
    torch::torch_manual_seed(13)

    net <- torch::nn_module(
      "cnn",
      initialize = function() {
        self$conv1 <- torch::nn_conv2d(in_channels = 7, out_channels = 14, kernel_size = 3, padding = 0)
        self$conv2 <- torch::nn_conv2d(in_channels = 14, out_channels = 28, kernel_size = 3, padding = 0)
        self$fc1 <- torch::nn_linear(in_features = 7 * 7 * 28, out_features = 28)
        self$fc2 <- torch::nn_linear(in_features = 28, out_features = 1)
      },
      forward = function(x) {
        x %>%
          self$conv1() %>%
          torch::nnf_relu() %>%
          self$conv2() %>%
          torch::nnf_relu() %>%
          torch::torch_flatten(start_dim = 2) %>%
          self$fc1() %>%
          torch::nnf_relu() %>%
          self$fc2()
      }
    )
    ##

    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      message("-- Evaluating with fold ", j, "/", length(folds))

      train_dataloader <-
        data[data[, partition] != folds[j], c(longitude, latitude, response)] %>%
        cnn_make_samples(longitude, latitude, response, rasters, crop_size) %>%
        create_dataset() %>%
        torch::dataloader(batch_size = batch_size, shuffle = TRUE)

      test_dataloader <-
        data[data[, partition] == folds[j], c(longitude, latitude, response)] %>%
        cnn_make_samples(longitude = longitude, latitude = latitude, response = response, raster = rasters, size = crop_size) %>%
        create_dataset() %>%
        torch::dataloader(batch_size = batch_size, shuffle = TRUE)

      # fit model
      model <- net %>%
        luz::setup(
          loss = torch::nn_l1_loss(),
          optimizer = optim_adam
        ) %>%
        luz::set_opt_hparams(lr = learning_rate) %>%
        fit(train_dataloader, epochs = n_epochs, valid_data = test_dataloader)

      pred <- predict(model, test_dataloader) %>% as.numeric()
      observed <- test_dataloader$dataset$response_variable %>% as.numeric()
      eval_partial[[j]] <- dplyr::tibble(
        model = "cnn",
        adm_eval(obs = observed, pred = pred)
      )

      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }

    # fit final model with all data

    # nota: precisa criar um torch dataset e um dataloader para todos os dados

    full_dataloader <- data[, c(longitude, latitude, response)] %>%
      cnn_make_samples(longitude, latitude, response, rasters, crop_size) %>%
      create_dataset() %>%
      torch::dataloader(batch_size = batch_size, shuffle = TRUE)

    # nota: na sequência abaixo ele fitta um modelo com todos os dados
    full_model <- net %>%
      luz::setup(
        loss = torch::nn_l1_loss(),
        optimizer = optim_adam
      ) %>%
      luz::set_opt_hparams(lr = learning_rate) %>%
      fit(full_dataloader, epochs = n_epochs)

    # bind predicted evaluation
    eval_partial <- eval_partial %>%
      dplyr::bind_rows() %>%
      dplyr::as_tibble()

    # bind predicted partition
    if (predict_part) {
      part_pred <- part_pred %>%
        dplyr::bind_rows() %>%
        dplyr::as_tibble()
    } else {
      part_pred <- NULL
    }

    # Summarize performance
    eval_final <- eval_partial %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(dplyr::across(corr_spear:pdispersion, list(
        mean = mean,
        sd = stats::sd
      )), .groups = "drop")

    # Final object
    data_list <- list(
      model = full_model,
      predictors = variables,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )

    return(data_list)
  }

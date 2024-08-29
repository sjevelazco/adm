#' Fit and validate Convolutional Neural Network Model
#'
#' @description This function is used to fit a convolutional neural network (CNN) model for abundance.
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param x character. The name of the column containing longitude information for each observation.
#' @param y character. The name of the column containing latitude information for each observation.
#' @param rasters a terra SpatRaster object. A raster containing the predictor variables to be cropped around each observation.
#' @param crop_size numeric. An integer defining the range of pixels around the observation from the raster object passed to rasters parameter. Default = 5
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param learning_rate numeric. The size of the step taken during the optimization process. Default = 0.01
#' @param n_epochs numeric. How many times the learning algorithm will work through the training set. Default = 10
#' @param batch_size numeric. A batch is a subset of the training set used in a single iteration of the training process. The size of each batch is referred to as the batch size. Default = 32
#' @param custom_architecture a Torch nn_module_generator object. A neural network architecture to be used instead of the internal default one. Default NULL
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom luz setup set_opt_hparams fit
#' @importFrom stats sd
#' @importFrom terra rast
#' @importFrom torch dataset torch_tensor torch_manual_seed nn_module nn_conv2d nn_linear nnf_relu torch_flatten dataloader nn_l1_loss optim_adam
#' @importFrom torchvision transform_to_tensor
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "luz_module_fitted" object from luz (torch framework). This object can be used to predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
fit_abund_cnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           x,
           y,
           rasters,
           crop_size = 5,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           learning_rate = 0.01,
           n_epochs = 10,
           batch_size = 32,
           custom_architecture = NULL) {
    self <- corr_spear <- pdisp <- NULL
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
        response <- torch::torch_tensor(self$response_variable[[index]])
        x <- torchvision::transform_to_tensor(self$predictors[[index]])
        list(x = x, y = response)
      },
      .length = function() {
        length(self$response_variable)
      }
    )

    # loading rasters
    if (class(rasters) %in% "character") {
      rasters <- terra::rast(rasters)
      rasters <- rasters[[c(predictors, predictors_f)]]
    } else {
      stop("Please, provide a path to the raster file.")
    }

    ##
    torch::torch_manual_seed(13)

    if (!is.null(custom_architecture)) {
      net <- custom_architecture
    } else {
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
    }
    ##
    samples_list <- list()
    for (fold in folds) {
      fold_mtx <- data[data[, partition] == fold, c(x, y, response)] %>%
        cnn_make_samples(x, y, response, rasters, crop_size) %>%
        list()

      names(fold_mtx) <- fold

      samples_list <- append(samples_list, fold_mtx)
    }

    ##
    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      message("-- Partition number ", j, "/", length(folds))

      train_samples <- samples_list[folds[folds != j]]
      train_samples <- unlist(train_samples, recursive = FALSE)

      train_response_folds <- names(train_samples)[grep("\\.response$", names(train_samples))]
      train_concatened_responses <- list()
      for (fold in train_response_folds) {
        train_concatened_responses <- c(train_concatened_responses, train_samples[[fold]])
      }

      train_predictors_folds <- names(train_samples)[grep("\\.predictors$", names(train_samples))]
      train_concatened_predictors <- list()
      for (fold in train_predictors_folds) {
        train_concatened_predictors <- c(train_concatened_predictors, train_samples[[fold]])
      }


      train_data_list <- list(
        predictors = train_concatened_predictors,
        response = train_concatened_responses
      )

      train_dataloader <- create_dataset(train_data_list) %>%
        torch::dataloader(batch_size = batch_size, shuffle = TRUE)


      test_samples <- samples_list[folds[folds == j]]
      test_samples <- unlist(test_samples, recursive = FALSE)

      test_response_folds <- names(test_samples)[grep("\\.response$", names(test_samples))]
      test_concatened_responses <- list()
      for (fold in test_response_folds) {
        test_concatened_responses <- c(test_concatened_responses, test_samples[[fold]])
      }

      test_predictors_folds <- names(test_samples)[grep("\\.predictors$", names(test_samples))]
      test_concatened_predictors <- list()
      for (fold in test_predictors_folds) {
        test_concatened_predictors <- c(test_concatened_predictors, test_samples[[fold]])
      }


      test_data_list <- list(
        predictors = test_concatened_predictors,
        response = test_concatened_responses
      )

      test_dataloader <- create_dataset(test_data_list) %>%
        torch::dataloader(batch_size = batch_size, shuffle = TRUE)

      # fit model
      suppressMessages(
        model <- net %>%
          luz::setup(
            loss = torch::nn_l1_loss(),
            optimizer = torch::optim_adam
          ) %>%
          luz::set_opt_hparams(lr = learning_rate) %>%
          luz::fit(train_dataloader, epochs = n_epochs, valid_data = test_dataloader)
      )

      pred <- predict(model, test_dataloader) # cuda atualization
      pred <- pred$to(device = "cpu") #
      pred <- as.numeric(pred) #
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


    full_samples <- unlist(samples_list, recursive = FALSE)

    full_response_folds <- names(full_samples)[grep("\\.response$", names(full_samples))]
    full_concatened_responses <- list()
    for (fold in full_response_folds) {
      full_concatened_responses <- c(full_concatened_responses, full_samples[[fold]])
    }

    full_predictors_folds <- names(full_samples)[grep("\\.predictors$", names(full_samples))]
    full_concatened_predictors <- list()
    for (fold in full_predictors_folds) {
      full_concatened_predictors <- c(full_concatened_predictors, full_samples[[fold]])
    }


    full_data_list <- list(
      predictors = full_concatened_predictors,
      response = full_concatened_responses
    )

    full_dataloader <- create_dataset(full_data_list) %>%
      torch::dataloader(batch_size = batch_size, shuffle = TRUE)

    # nota: na sequência abaixo ele fitta um modelo com todos os dados
    suppressMessages(
      full_model <- net %>%
        luz::setup(
          loss = torch::nn_l1_loss(),
          optimizer = torch::optim_adam
        ) %>%
        luz::set_opt_hparams(lr = learning_rate) %>%
        luz::fit(full_dataloader, epochs = n_epochs)
    )

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
      dplyr::summarise(dplyr::across(corr_spear:pdisp, list(
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

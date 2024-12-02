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
#' @param sample_size numeric. A vector containing the dimensions, in pixels, of raster samples. See cnn_make_samples beforehand. Default c(11,11)
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param learning_rate numeric. The size of the step taken during the optimization process. Default = 0.01
#' @param n_epochs numeric. Maximum number of times the learning algorithm will work through the training set. Default = 10
#' @param batch_size numeric. A batch is a subset of the training set used in a single iteration of the training process. The size of each batch is referred to as the batch size. Default = 32
#' @param custom_architecture a Torch nn_module_generator object. A neural network architecture to be used instead of the internal default one. Default NULL
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#' @param validation_patience numerical. An integer indicating the number of epochs without loss improvement tolerated by the algorithm in the validation process. If the patience limit is exceeded, the training ends. Default 2
#' @param fitting_patience numerical. The same as validation_patience, but in the final model fitting process. Default 5
#'
#' @importFrom dplyr bind_rows select starts_with pull tibble as_tibble group_by summarise across bind_cols
#' @importFrom luz setup set_opt_hparams fit luz_callback_early_stopping
#' @importFrom stats sd
#' @importFrom terra rast
#' @importFrom torch dataset torch_tensor torch_manual_seed dataloader nn_l1_loss optim_adam
#' @importFrom torchvision transform_to_tensor
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "luz_module_fitted" object from luz (torch framework). This object can be used to predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each replica and partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Database with species abundance and x and y coordinates
#' data("sppabund")
#'
#' # Extract data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(-.part2, -.part3)
#'
#' # Explore reponse variables
#' some_sp$ind_ha %>% range()
#' some_sp$ind_ha %>% hist()
#'
#' # Here we balance number of absences
#' some_sp <-
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
#'
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)
#'
#' # Generate an architecture
#' cnn_arch <- generate_cnn_architecture(
#'   number_of_features = 3,
#'   number_of_outputs = 1,
#'   sample_size = c(11, 11),
#'   number_of_conv_layers = 2,
#'   conv_layers_size = c(14, 28),
#'   conv_layers_kernel = 3,
#'   conv_layers_stride = 1,
#'   conv_layers_padding = 0,
#'   number_of_fc_layers = 1,
#'   fc_layers_size = c(28),
#'   pooling = NULL,
#'   batch_norm = TRUE,
#'   dropout = 0,
#'   verbose = T
#' )
#'
#' # Fit a CNN model
#' mcnn <- fit_abund_cnn(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = NULL,
#'   partition = ".part",
#'   x = "x",
#'   y = "y",
#'   rasters = envar,
#'   sample_size = c(11, 11),
#'   learning_rate = 0.01,
#'   n_epochs = 100,
#'   batch_size = 32,
#'   validation_patience = 2,
#'   fitting_patience = 5,
#'   custom_architecture = cnn_arch,
#'   verbose = TRUE,
#'   predict_part = TRUE
#' )
#'
#' mcnn
#' }
fit_abund_cnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           x,
           y,
           rasters,
           sample_size,
           partition,
           predict_part = FALSE,
           learning_rate = 0.01,
           n_epochs = 10,
           batch_size = 32,
           validation_patience = 2,
           fitting_patience = 5,
           custom_architecture = NULL,
           verbose = TRUE) {
    # TODO check argument description validation_patience and fitting_patience
    . <- self <- corr_spear <- pdisp <- envar <- mae <- NULL
    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- dplyr::bind_rows(c(c = predictors))
    }

    # Adequate database
    data <- adapt_df(
      data = data,
      predictors = predictors,
      predictors_f = predictors_f,
      response = response,
      partition = partition,
      xy = c("x", "y")
    )

    # Get cropsize
    crop_size <- cnn_get_crop_size(sample_size = sample_size)

    # # ---- Formula ----
    # if (is.null(fit_formula)) {
    #   formula1 <- stats::formula(paste(response, "~", paste(c(
    #     predictors,
    #     predictors_f
    #   ), collapse = " + ")))
    # } else {
    #   formula1 <- fit_formula
    # }

    # if (verbose) {
    #   message(
    #     "Formula used for model fitting:\n",
    #     Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
    #     "\n"
    #   )
    # }


    # create_dataset definition
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
    } else if (class(rasters) %in% "SpatRaster") {
      rasters <- rasters[[c(predictors, predictors_f)]]
    } else {
      stop("Please, provide a SpatRaster object or a path to the raster file.")
    }

    # architecture setup
    torch::torch_manual_seed(13)

    if (!is.null(custom_architecture)) {
      if ("arch" %in% names(custom_architecture)) {
        custom_architecture <- custom_architecture$net
      }
      net <- custom_architecture
    } else {
      net <- generate_cnn_architecture(
        number_of_features = length(variables),
        number_of_outputs = 1,
        sample_size = rep((crop_size * 2) + 1, 2),
        number_of_conv_layers = 2,
        conv_layers_size = c(length(variables), length(variables)),
        conv_layers_kernel = 3,
        conv_layers_stride = 1,
        conv_layers_padding = 0,
        number_of_fc_layers = 1,
        fc_layers_size = c(length(variables)),
        pooling = NULL,
        batch_norm = TRUE,
        dropout = 0,
        verbose = FALSE
      )$net
    }

    # Fit models
    np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))

    part_pred_list <- list()
    eval_partial_list <- list()

    for (h in 1:np) {
      if (verbose) {
        message("Replica number: ", h, "/", np)
      }
      # out <- pre_tr_te(data, p_names, h)

      folds <- data %>%
        dplyr::pull(p_names[h]) %>%
        unique() %>%
        sort()

      samples_list <- list()
      for (fold in folds) {
        fold_mtx <- data[data[, p_names[h]] == fold, c(x, y, response)] %>%
          cnn_make_samples(x, y, response, rasters, size = crop_size) %>%
          list()

        names(fold_mtx) <- fold

        samples_list <- append(samples_list, fold_mtx)
      }
      rm(fold_mtx)

      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()

      for (j in 1:length(folds)) {
        if (verbose) {
          message("-- Partition number ", j, "/", length(folds))
        }

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
        rm(
          train_concatened_predictors,
          train_concatened_responses
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
        set.seed(13)
        # suppressMessages(
        model <- net %>%
          luz::setup(
            loss = torch::nn_l1_loss(),
            optimizer = torch::optim_adam
          ) %>%
          luz::set_opt_hparams(lr = learning_rate) %>%
          luz::fit(train_dataloader,
            epochs = n_epochs,
            valid_data = test_dataloader,
            callbacks = luz::luz_callback_early_stopping(patience = validation_patience)
          )
        # )

        pred <- predict(model, test_dataloader)
        pred <- pred$to(device = "cpu")
        pred <- as.numeric(pred)
        observed <- test_dataloader$dataset$response_variable %>% as.numeric()
        eval_partial[[j]] <- dplyr::tibble(
          model = "cnn",
          adm_eval(obs = observed, pred = pred)
        )

        if (predict_part) {
          part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
        }
      }

      # Create final database with parameter performance
      names(eval_partial) <- 1:length(folds)
      eval_partial <-
        eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial

      if (predict_part) {
        names(part_pred) <- 1:length(folds)
        part_pred <-
          part_pred[sapply(part_pred, function(x) !is.null(dim(x)))] %>%
          dplyr::bind_rows(., .id = "partition")
        part_pred_list[[h]] <- part_pred
      }
    }

    # fit final model with all data

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

    set.seed(13)
    # suppressMessages(
    full_model <- net %>%
      luz::setup(
        loss = torch::nn_l1_loss(),
        optimizer = torch::optim_adam
      ) %>%
      luz::set_opt_hparams(lr = learning_rate) %>%
      luz::fit(full_dataloader,
        epochs = n_epochs,
        callbacks = luz::luz_callback_early_stopping(
          monitor = "train_loss",
          patience = fitting_patience
        )
      )
    # )

    # bind predicted evaluation
    eval_partial <- eval_partial_list %>%
      dplyr::bind_rows(.id = "replica") %>%
      dplyr::as_tibble()

    # bind predicted partition
    if (predict_part) {
      part_pred <- part_pred_list %>%
        dplyr::bind_rows(.id = "replica")
    } else {
      part_pred <- NULL
    }

    # Sumarize performance
    eval_final <- eval_partial %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(
        dplyr::across(
          c(mae:pdisp),
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ stats::sd(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      )

    variables <- dplyr::bind_cols(
      data.frame(
        model = "cnn",
        response = response
      ),
      variables
    ) %>% as_tibble()

    # Final object
    data_list <- list(
      model = full_model,
      predictors = variables,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )

    # Standardize output list
    for (i in 2:length(data_list)) {
      if (!class(data_list[[i]])[1] == "tbl_df") {
        data_list[[i]] <- dplyr::as_tibble(data_list[[i]])
      }
    }

    return(data_list)
  }

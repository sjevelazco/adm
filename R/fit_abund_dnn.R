#' Fit and validate Deep Neural Network model
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param learning_rate numeric. The size of the step taken during the optimization process. Default = 0.01
#' @param n_epochs numeric. Max number of times the learning algorithm will work through the training set. Default = 10
#' @param batch_size numeric. A batch is a subset of the training set used in a single iteration of the training process. The size of each batch is referred to as the batch size. Default = 32
#' @param custom_architecture a Torch nn_module_generator object or a generate_dnn_architecture output. A neural network architecture to be used instead of the internal default one. Default NULL
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#' @param validation_patience numerical. An integer indicating the number of epochs without loss improvement tolerated by the algorithm in the validation process. If the patience limit is exceeded, the training ends. Default 2
#' @param fitting_patience numerical. The same as validation_patience, but in the final model fitting process. Default 5
#'
#' @importFrom dplyr bind_rows bind_cols pull tibble as_tibble group_by summarise across
#' @importFrom luz setup set_opt_hparams fit
#' @importFrom stats sd runif
#' @importFrom torch dataset torch_tensor torch_manual_seed nn_module nn_linear nnf_relu dataloader nn_l1_loss optim_adam
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
#' # Generate a architecture
#' dnn_arch <- generate_dnn_architecture(
#'   number_of_features = 3,
#'   number_of_outputs = 1,
#'   number_of_hidden_layers = 3,
#'   hidden_layers_size = c(8, 16, 8),
#'   batch_norm = TRUE
#' )
#'
#' # Fit a NET model
#' mdnn <- fit_abund_dnn(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = NULL,
#'   partition = ".part",
#'   learning_rate = 0.01,
#'   n_epochs = 10,
#'   batch_size = 32,
#'   validation_patience = 2,
#'   fitting_patience = 5,
#'   custom_architecture = dnn_arch,
#'   verbose = TRUE,
#'   predict_part = TRUE
#' )
#'
#' mdnn
#' }
fit_abund_dnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
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

    . <- self <- model <- mae <- pdisp <- NULL
    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- dplyr::bind_rows(c(c = predictors))
    }

    if (!is.null(predictors_f)) {
      warning("Categorical variables aren't available for DNN and will be ignored.")
      predictors_f <- NULL
    }

    # Adequate database
    data <- adapt_df(
      data = data,
      predictors = predictors,
      predictors_f = predictors_f,
      response = response,
      partition = partition
    )

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
      initialize = function(df, response_variable) {
        self$df <- df[, -which(names(df) == response_variable)]
        self$response_variable <- df[[response_variable]]
      },
      .getitem = function(index) {
        response <- torch::torch_tensor(self$response_variable[index])
        x <- torch::torch_tensor(as.numeric(self$df[index, ]))
        list(x = x, y = response)
      },
      .length = function() {
        length(self$response_variable)
      }
    )

    # architecture setup
    torch::torch_manual_seed(13)

    if (!is.null(custom_architecture)) {
      if ("arch" %in% names(custom_architecture)) {
        custom_architecture <- custom_architecture$net
      }
      net <- custom_architecture
    } else {
      net <- torch::nn_module(
        "neural_net",
        initialize = function() {
          self$input <- torch::nn_linear(length(variables), 2 * length(variables))
          self$linear1 <- torch::nn_linear(2 * length(variables), length(variables))
          self$output <- torch::nn_linear(length(variables), 1)
        },
        forward = function(x) {
          x %>%
            self$input() %>%
            torch::nnf_relu() %>%
            self$linear1() %>%
            torch::nnf_relu() %>%
            self$output()
        }
      )
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

      folds <- data %>%
        dplyr::pull(p_names[h]) %>%
        unique() %>%
        sort()

      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()

      for (j in 1:length(folds)) {
        if (verbose) {
          message("-- Partition number ", j, "/", length(folds))
        }

        train_set <- data[data[, p_names[h]] != folds[j], c(predictors, response)] %>%
          create_dataset(response_variable = response)
        test_set <- data[data[, p_names[h]] == folds[j], c(predictors, response)] %>%
          create_dataset(response_variable = response)

        train_dataloader <- torch::dataloader(train_set, batch_size = batch_size, shuffle = TRUE)
        test_dataloader <- torch::dataloader(test_set, batch_size = batch_size, shuffle = TRUE)

        set.seed(13)
        suppressMessages(
          fitted <- net %>%
            luz::setup(
              loss = torch::nn_l1_loss(),
              optimizer = torch::optim_adam
            ) %>%
            luz::set_opt_hparams(lr = learning_rate) %>%
            luz::fit(train_dataloader,
                     valid_data = test_dataloader,
                     epochs = n_epochs,
                     callbacks = luz::luz_callback_early_stopping(patience = validation_patience)
            )
        )
        
        pred <- predict(fitted, test_set) %>% as.numeric()

        if (!(sum(is.na(pred)) == length(pred))) {
          pred[is.na(pred)] <- stats::runif(sum(is.na(pred)), min(data[[response]]), max(data[[response]]))
          observed <- test_set$response_variable %>% as.numeric()
          eval_partial[[j]] <- dplyr::tibble(
            model = "dnn",
            adm_eval(obs = observed, pred = pred)
          )
        }

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
    df <- create_dataset(data[, c(predictors, response)], response)
    df_dl <- torch::dataloader(df, batch_size = batch_size, shuffle = TRUE)

    set.seed(13)
    suppressMessages(
      full_fitted <- net %>%
        luz::setup(
          loss = torch::nn_l1_loss(),
          optimizer = torch::optim_adam
        ) %>%
        luz::set_opt_hparams(lr = learning_rate) %>%
        luz::fit(df_dl,
                 epochs = n_epochs,
                 callbacks = luz::luz_callback_early_stopping(
                   monitor = "train_loss",
                   patience = fitting_patience
                 )
        )
    )
    
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
            sd = ~ sd(.x, na.rm = TRUE)
          )
        ),
        .groups = "drop"
      )

    variables <- dplyr::bind_cols(
      data.frame(
        model = "dnn",
        response = response
      ),
      variables
    ) %>% as_tibble()

    # Final object
    data_list <- list(
      model = full_fitted,
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

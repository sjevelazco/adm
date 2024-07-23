#' Fit and validate Deep Neural Network model
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
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
#' @importFrom torch dataset torch_tensor torch_manual_seed nn_module nn_linear nnf_relu dataloader nn_l1_loss optim_adam
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
fit_abund_dnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           learning_rate = 0.01,
           n_epochs = 10,
           batch_size = 32,
           custom_architecture = NULL) {
    # Variables
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()

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

    ##
    torch::torch_manual_seed(13)

    if (!is.null(custom_architecture)) {
      net <- custom_architecture
    } else {
      net <- torch::nn_module(
        "neural_net",
        initialize = function() {
          self$input <- torch::nn_linear(ncol(variables), 2 * ncol(variables))
          self$linear1 <- torch::nn_linear(2 * ncol(variables), ncol(variables))
          self$output <- torch::nn_linear(ncol(variables), 1)
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
    ##

    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      message("-- Evaluating with fold ", j, "/", length(folds))

      # nota: nesta parte se cria dois torch datasets, um para treino, outro para teste
      train_set <- data[data[, partition] != folds[j], c(predictors, response)] %>%
        create_dataset(response_variable = response)
      test_set <- data[data[, partition] == folds[j], c(predictors, response)] %>%
        create_dataset(response_variable = response)

      train_dataloader <- torch::dataloader(train_set, batch_size = batch_size, shuffle = TRUE)
      test_dataloader <- torch::dataloader(test_set, batch_size = batch_size, shuffle = TRUE)

      fitted <- net %>%
        luz::setup(
          loss = torch::nn_l1_loss(),
          optimizer = torch::optim_adam
        ) %>%
        luz::set_opt_hparams(lr = learning_rate) %>%
        luz::fit(train_dataloader, epochs = n_epochs)
      # valid_data = test_dataloader)

      pred <- predict(fitted, test_set) %>% as.numeric() # nota: não existe mais objeto "model", agora é "fitted"

      if (!(sum(is.na(pred)) == length(pred))) {
        pred[is.na(pred)] <- runif(sum(is.na(pred)), 0, 100)
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

    # fit final model with all data

    df <- create_dataset(data[, c(predictors, response)], response)

    df_dl <- torch::dataloader(df, batch_size = batch_size, shuffle = TRUE)

    # nota: na sequência abaixo ele fitta um modelo com todos os dados
    full_fitted <- net %>%
      luz::setup(
        loss = torch::nn_l1_loss(),
        optimizer = torch::optim_adam
      ) %>%
      luz::set_opt_hparams(lr = learning_rate) %>%
      luz::fit(df_dl, epochs = n_epochs)

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
      model = full_fitted,
      predictors = variables,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )

    return(data_list)
  }

#' Fit and validate Artificial Neural Network models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param size numerical. The size of the hidden layer.
#' @param decay numerial. Value for weight decay. Default 0.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom dplyr pull tibble bind_rows as_tibble group_by summarise across
#' @importFrom kernlab ksvm predict
#' @importFrom stats formula sd
#' @importFrom nnet nnet
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "ksvm" class object from kernlab package. This object can be used for predicting.
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
#' # Fit a NET model
#' mnet <- fit_abund_net(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   size = 32,
#'   decay = 0.1,
#'   predict_part = TRUE
#' )
#'
#' mnet
#' }
fit_abund_net <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           size,
           decay = 0,
           verbose = TRUE) {
    . <- mae <- pdisp <- NULL

    # Adequate database
    data <- adapt_df(
      data = data,
      predictors = predictors,
      predictors_f = predictors_f,
      response = response,
      partition = partition
    )

    # Variables
    variables <- get_variables(predictors, predictors_f)

    # Formula
    formula1 <- infer_formula(fit_formula, response, predictors, predictors_f, verbose)

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
        train_set <- data[data[, p_names[h]] != folds[j], ]
        test_set <- data[data[, p_names[h]] == folds[j], ]

        set.seed(13)
        model <- nnet::nnet(
          formula1,
          data = train_set,
          size = size,
          decay = decay,
          rang = 0.1,
          maxit = 1000,
          reltol = 1e-5,
          linout = TRUE,
          trace = FALSE
        )

        pred <- suppressMessages(stats::predict(model, newdata = test_set, type = "raw"))
        observed <- dplyr::pull(test_set, response)
        eval_partial[[j]] <- dplyr::tibble(
          model = "net",
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
    set.seed(13)
    full_model <- nnet::nnet(
      formula1,
      data = data,
      size = size,
      decay = decay,
      rang = 0.1,
      maxit = 1000,
      reltol = 1e-5,
      linout = TRUE,
      trace = FALSE
    )

    # Construct the standard final list to be returned
    data_list <- wrap_final_list(
      "net",
      full_model, 
      variables, 
      response, 
      eval_partial_list, 
      predict_part, 
      part_pred_list,
      get_metadata(
        "net", 
        list(
          formula = formula1,
          rang = 0.1,
          maxit = 1000,
          reltol = 1e-5,
          linout = "TRUE",
          trace = "FALSE"
        )
      )
    )

    return(data_list)
  }

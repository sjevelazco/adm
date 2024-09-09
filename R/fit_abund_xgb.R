#' Fit and validate Extreme Gradient Boosting models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE.
#' @param nrounds integer. Max number of boosting iterations. Default is 100.
#' @param max_depth integer. The maximum depth of each tree. Default 5
#' @param eta numeric. The learning rate of the algorithm. Default 0.1
#' @param gamma numeric. Minimum loss reduction required to make a further partition on a leaf node of the tree. Default is 1.
#' @param colsample_bytree numeric. Subsample ratio of columns when constructing each tree. Default is 1.
#' @param min_child_weight numeric. Minimum sum of instance weight needed in a child. Default is 1.
#' @param subsample numeric. Subsample ratio of the training instance. Default is 0.5.
#' @param objective character. The learning task and the corresponding learning objective. Default is "reg:squarederror", regression with squared loss.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE.
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom stats predict sd
#' @importFrom xgboost xgboost
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "xgb.Booster" class object from xgboost package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each replica and partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
fit_abund_xgb <-
  function(data,
           response,
           predictors,
           predictors_f,
           partition,
           predict_part = FALSE,
           nrounds = 100,
           max_depth = 5,
           eta = 0.1,
           gamma = 1,
           colsample_bytree = 1,
           min_child_weight = 1,
           subsample = 0.5,
           objective = "reg:squarederror",
           verbose = TRUE) {
    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- predictors
    }
    
    # Adequate database
    data <- adapt_df(data = data,
                     predictors = predictors,
                     predictors_f = predictors_f,
                     response = response,
                     partition = partition)

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
      
      folds <- data %>% dplyr::pull(p_names[h]) %>% unique() %>% sort()
      
      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()
      
      for (j in 1:length(folds)) {
        if (verbose) {
          message("-- Partition number ", j, "/", length(folds))
        }
        train_set <- data[data[, p_names[h]] != folds[j], ]
        test_set <- data[data[, p_names[h]] == folds[j], ]
        
        
        sp_train <- list(
          data = as.matrix(train_set[, variables]),
          target = train_set[, response]
        )
        
        sp_test <- list(
          data = as.matrix(test_set[, variables]),
          target = test_set[, response]
        )
        
        model <- xgboost::xgboost(
          data = sp_train$data,
          label = sp_train$target,
          params = list(
            max_depth = max_depth,
            eta = eta,
            gamma = gamma,
            colsample_bytree = colsample_bytree,
            min_child_weight = min_child_weight,
            subsample = subsample,
            objective = objective
          ),
          nrounds = nrounds,
          verbose = 0
        )
        
        pred <- stats::predict(model, sp_test$data, type = "response")
        observed <- sp_test$target
        eval_partial[[j]] <- dplyr::tibble(
          model = "xgb",
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
    full_model <- xgboost::xgboost(
      data = as.matrix(data[, variables]),
      label = data[, response],
      params = list(
        max_depth = max_depth,
        eta = eta,
        gamma = gamma,
        colsample_bytree = colsample_bytree,
        min_child_weight = min_child_weight,
        subsample = subsample,
        objective = objective
      ),
      nrounds = nrounds,
      verbose = 0
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

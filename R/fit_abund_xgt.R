fit_abund_xgt <-
  function(data,
           response,
           predictors,
           predictors_f,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           params,
           nrounds = 500,
           verbose = TRUE) {
    # Variables
    if (!is.null(predictors_f)){
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- predictors
    }
    
    # # ---- Formula ----
    # if (is.null(fit_formula)) {
    #   formula1 <- stats::formula(paste(response, "~", paste(c(
    #     predictors,
    #     predictors_f
    #   ), collapse = " + ")))
    # } else {
    #   formula1 <- fit_formula
    # }
    
    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()
    
    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      if(verbose){
        message("-- Evaluating with fold ", j, "/", length(folds))
      }
      
      train_set <- data[data[, partition] != folds[j], ]
      test_set <- data[data[, partition] == folds[j], ]
      
      sp_train <- list(
        data = as.matrix(train_set[, variables]),
        target = train_set[, response]
      )
      
      sp_test <- list(
        data = as.matrix(test_set[, variables]),
        target = test_set[, response]
      )
      
      #
      part_model <- xgboost::xgboost(
        data = sp_train$data,
        label = sp_train$target[[1]],
        params = params,
        nrounds = nrounds,
        verbose = verbose
      )
      #
      
      pred <- stats::predict(part_model, sp_test$data, type = "response")
      observed <- sp_test$target[[1]]
      eval_partial[[j]] <- dplyr::tibble(
        model = "xgt",
        adm_eval(obs = observed, pred = pred)
      )
      
      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }
    
    # fit final model with all data
    model <- xgboost::xgboost(
      data = as.matrix(data[, variables]),
      label = data[, response][[1]],
      params = params,
      nrounds = nrounds,
      verbose = verbose
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
      model = model,
      predictors = variables,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )
    return(data_list)
  }
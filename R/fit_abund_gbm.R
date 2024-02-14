#' Fit and validate Random Forests models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor
#' variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param mtry numeric. Number of variables randomly sampled as candidates at each split. Default sqrt(length(c(predictors, predictors_f)))
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "randomForest" class object from randomForest package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
fit_abund_gbm <-
  function(data,
           response,
           predictors,
           predictors_f,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           distribution,
           n.trees = 100,
           interaction.depth = 1,
           n.minobsinnode = 5,
           shrinkage = 0.1,
           verbose = TRUE) {
    # Variables
    if (!is.null(predictors_f)){
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- predictors
    }
    
    
    # ---- Formula ----
    if (is.null(fit_formula)) {
      formula1 <- stats::formula(paste(response, "~", paste(c(
        predictors,
        predictors_f
      ), collapse = " + ")))
    } else {
      formula1 <- fit_formula
    }
    
    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()
    
    # ---- Distribution ----
    if(distribution == "poisson"){
      data[,response] <- round(data[,response])
    }
    
    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      if(verbose){
        message("-- Evaluating with fold ", j, "/", length(folds))
      }
      
      train_set <- data[data[, partition] != folds[j], ]
      test_set <- data[data[, partition] == folds[j], ]
      
      #
      part_model <- gbm::gbm(
        formula = formula1,
        data = train_set,
        distribution = distribution,
        n.trees = n.trees,
        interaction.depth	= interaction.depth,
        n.minobsinnode = n.minobsinnode,
        shrinkage = shrinkage,
        bag.fraction = 0.9
      )
      
      #
      
      pred <- stats::predict(part_model, test_set, type = "response")
      observed <- test_set$ind_ha
      eval_partial[[j]] <- dplyr::tibble(
        model = "gbm",
        adm_eval(obs = observed, pred = pred)
      )
      
      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }
    
    # fit final model with all data
    model <- gbm::gbm(
      formula = formula1,
      data = data,
      distribution = distribution,
      n.trees = n.trees,
      interaction.depth	= interaction.depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = shrinkage,
      bag.fraction = 0.9
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

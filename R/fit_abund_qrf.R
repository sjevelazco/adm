#' Fit and validate Quantile Regression Random Forests models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param framework 
#' @param train_quantiles 
#' @param eval_quantile 
#' @param mtry numeric. Number of variables randomly sampled as candidates at each split. Default (length(c(predictors, predictors_f))/3)
#' @param ntree numeric. Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. Default 500
#' @param nodesize 
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom randomForest randomForest
#' @importFrom stats formula sd
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "randomForest" class object from randomForest package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each replica and partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#' @seealso \code{\link{fit_abund_cnn}},
#' \code{\link{fit_abund_dnn}},
#' \code{\link{fit_abund_gam}},
#' \code{\link{fit_abund_gbm}},
#' \code{\link{fit_abund_glm}},
#' \code{\link{fit_abund_net}},
#' \code{\link{fit_abund_svm}},
#' \code{\link{fit_abund_xgb}}
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
#' }
fit_abund_qrf <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           framework = "quantregForest",
           train_quantiles = c(0.5),
           eval_quantile = 0.5,
           mtry = length(c(predictors, predictors_f)) / 3,
           ntree = 2000,
           nodesize = 5,
           verbose = TRUE) {
    . <- mae <- pdisp <- NULL
    
    # Algorithm checks and fixes
    if(framework=="grf"){
      predictors_f <- NULL
    }
    
    if (!eval_quantile %in% train_quantiles){
      train_quantiles <- sort(c(train_quantiles, eval_quantile))
    }
    
    # Adequate database
    data <- adapt_df(data, predictors, predictors_f, response, partition)
    
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
        
        model <- switch (framework,
          "grf" = {
            grf::quantile_forest(
              mtry = mtry,
              num.trees = ntree,
              min.node.size = nodesize,
              X = train_set[,c(predictors)],
              Y = train_set[,response],
              quantiles = train_quantiles,
              num.threads = 1
            )
          },
          "quantregForest" = {
            quantregForest::quantregForest(
              formula = formula1,
              mtry = mtry,
              ntree = ntree,
              nodesize = nodesize,
              importance = FALSE,
              x = train_set[,c(predictors,predictors_f)],
              y = train_set[,response]
            )
          }
        )
        
        pred <- suppressMessages(
          stats::predict(
            model, 
            test_set[,c(predictors, predictors_f)], 
            type = "response", 
            what = train_quantiles
          )
        )
        
        pred <- switch (framework,
          "grf" = {
            pred$predictions[,which(train_quantiles == eval_quantile)]
          },
          "quantregForest" = {
            if (length(train_quantiles) > 1) {
              pred[,which(train_quantiles == eval_quantile)]
            } else {
              pred
            }
          }
        )
        
        observed <- dplyr::pull(test_set, response)
        
        eval_partial[[j]] <- dplyr::tibble(
          model = "qrf",
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
    full_model <- switch(framework,
      "grf" = {
        grf::quantile_forest(
          mtry = mtry,
          num.trees = ntree,
          X = data[, c(predictors)],
          Y = data[, response],
          quantiles = train_quantiles,
          num.threads = 1
        )
      },
      "quantregForest" = {
        quantregForest::quantregForest(
          formula = formula1,
          mtry = mtry,
          ntree = ntree,
          importance = FALSE,
          x = data[, c(predictors, predictors_f)],
          y = data[, response]
        )
      }
    )
    
    # Construct the standard final list to be returned
    data_list <- wrap_final_list(
      "qrf",
      full_model, 
      variables, 
      response, 
      eval_partial_list, 
      predict_part, 
      part_pred_list,
      get_metadata(
        "qrf", 
        list(
          train_quantiles = train_quantiles, 
          eval_quantile = eval_quantile, 
          framework = framework
        )
      )
    )
    
    return(data_list)
  }

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
#' 
#' @importFrom dplyr pull tibble bind_rows as_tibble group_by summarise across
#' @importFrom kernlab ksvm predict
#' @importFrom stats formula sd
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
    
    
    # Adequate database
    data <- adapt_df(data = data,
                     predictors = predictors,
                     predictors_f = predictors_f,
                     response = response,
                     partition = partition)
    
    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- predictors
    }
    
    # Formula
    if (is.null(fit_formula)) {
      formula1 <- stats::formula(paste(response, "~", paste(c(
        predictors,
        predictors_f
      ), collapse = " + ")))
    } else {
      formula1 <- fit_formula
    }
    
    if (verbose) {
      message(
        "Formula used for model fitting:\n",
        Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
        "\n"
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
      dplyr::summarise(dplyr::across(mae:pdisp, list(
        mean = mean,
        sd = stats::sd
      )), .groups = "drop")
    
    # Final object
    data_list <- list(
      model = full_model,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )
    return(data_list)
  }

#' fit_abund_gam
#'
#' @param data 
#' @param response 
#' @param predictors 
#' @param predictors_f 
#' @param fit_formula 
#' @param partition 
#' @param predict_part 
#' @param k 
#' @param family 
#' @param method 
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom mgcv gam
#' @importFrom stats formula predict sd
#'
#' @return
#' @export
#'
#' @examples
fit_abund_gam <-
  function(data,
           response,
           predictors,
           predictors_f,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           k = 10,
           family = "poisson",
           method = "REML") {
    
    # Variables
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    
    # Formula
    if (is.null(fit_formula)) {
      formula1 <-
        paste(c(
          paste("s(", predictors, paste0(", k = ", k, ")"), collapse = " + ", sep = ""),
          predictors_f
        ), collapse = " + ")
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    } else {
      formula1 <- fit_formula
    }
    
    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()
    
    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      message("-- Evaluating with fold ", j, "/", length(folds))
      
      data[,response] <- as.integer(round(data[,response],0))
      
      train_set <- data[data[, partition] != folds[j], ] 
      test_set <- data[data[, partition] == folds[j], ]
      
      model <- mgcv::gam(formula = formula1,
                         family = family,
                         data = train_set,
                         method = method)
      
      pred <- stats::predict(model, test_set, type = "response")
      observed <- dplyr::pull(test_set, response)
      eval_partial[[j]] <- dplyr::tibble(
        model = "gam",
        adm_eval(obs = observed, pred = pred)
      )
      
      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }
    
    # fit final model with all data
    full_model <- mgcv::gam(formula = formula1,
                            family = family,
                            data = train_set,
                            method = method)
    
    
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
      dplyr::summarise(dplyr::across(corr_spear:pdispersion, list(
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

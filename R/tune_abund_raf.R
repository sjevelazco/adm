#' Fit and validate Random Forest models with exploration of hyper-parameters that optimize performance
#'
#' @param data data.frame. Database with response (i.e., abundance) and predictors values.
#' @param response character. Column name with species abundance data (e.g., 0, 1, 45).
#' @param predictors character. Vector with the column names of quantitative predictor variables
#' (i.e. continuous variables). Usage predictors = c("temperature", "sand", "elevation")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables
#' (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(abund ~ temp + sand + pH + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default NULL
#' @param partition  character. Column name with training and validation partition groups.
#' @param predict_part 
#' @param grid  data.frame. A data frame object with algorithm hyper-parameters values to be tested.
#' It is recommended to generate this data.frame with the grid() function. Hyper-parameter needed
#' for tuning is 'mtry'. The maximum mtry cannot exceed the total number of predictors.
#' @param metrics character. Performance metric used for selecting the best combination of hyper-parameter values. One of the following metrics can be used: xxx, xxx, xxx, xxx, xxx,
#' AUC, and BOYCE. TSS is used as default.
#' @param n_cores numeric. Number of cores use for parallelization. Default 1
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
tune_abund_raf <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           metrics,
           n_cores = 1,
           verbose = FALSE) {
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Random Forest.")
      mtry <- seq(from = 2, to = length(predictors), by = 1)
      ntree <- seq(from = 500, to = 1000, by = 100)
      grid <- expand.grid(mtry = mtry, ntree = ntree)
    } else {
      if (all(names(grid)%in%c("mtry","ntree")) & length(names(grid))==2){
        grid <- grid
      } else {
        stop("Grid names expected to be mtry and ntree.")
      }
    }
    
    grid$comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    hyper_combinations <- list()
    for (i in 1:nrow(grid)) {
      cat(i,"/",nrow(grid)) # DEBUG
      model <-
        fit_abund_raf(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          mtry = grid[i,"mtry"],
          ntree = grid[i,"ntree"]
        )
      hyper_combinations <- append(hyper_combinations,list(model$performance))
      names(hyper_combinations)[i] <- grid[i,"comb_id"]
    }
    
    hyper_combinations <- dplyr::bind_rows(hyper_combinations, .id = "comb_id")

    ranked_combinations <- model_selection(hyper_combinations,metrics) %>%
      dplyr::left_join(grid, by = c("comb_id"="comb_id"))
    
    for (i in 1:length(ranked_combinations)) {
      ranked_combinations[[i]] <- ranked_combinations[[i]] %>%
        dplyr::left_join(grid, by = c("comb_id"="comb_id"))
    }
    
    # fit final model
    message("Fitting the best model...")
    final_model <-
      fit_abund_raf(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        mtry = ranked_combinations[[1]][[1,"mtry"]],
        ntree = ranked_combinations[[1]][[1,"ntree"]]
      )
    
    message('The best model was a Random Forest with mtry = ',
            ranked_combinations[[1]][[1,"mtry"]],
            " and ntree = ", 
            ranked_combinations[[1]][[1,"ntree"]])
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }

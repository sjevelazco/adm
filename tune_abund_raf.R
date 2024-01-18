#
tune_abund_raf <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           metrics) {
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Random Forest.")
      mtry <- seq(from = 1, to = length(predictors), by = 1)
      ntree <- seq(from = 50, to = 1000, by = 50)
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
      model <-
        fit_abund_raf(
          data = data,
          response = response,
          predictors = pred_var,
          predictors_f = NULL,
          fit_formula = NULL,
          partition = partition,
          predict_part = FALSE,
          mtry = grid[i,"mtry"],
          ntree = grid[i,"ntree"],
          custom_name = grid[i,"comb_id"]
        )
      hyper_combinations <- append(hyper_combinations,list(model$performance))
    }
    
    hyper_combinations <- do.call(rbind,hyper_combinations)

    ranked_combinations <- model_selection(hyper_combinations,metrics) %>%
      dplyr::left_join(grid, by = c("model"="comb_id"))
    
    # fit final model
    message("Fitting the best model...")
    final_model <-
      fit_abund_raf(
        data = data,
        response = response,
        predictors = pred_var,
        predictors_f = NULL,
        fit_formula = NULL,
        partition = partition,
        predict_part = FALSE,
        mtry = ranked_combinations[[1,"mtry"]],
        ntree = ranked_combinations[[1,"ntree"]],
        custom_name = "raf"
      )
    
    message('The best model was a Random Forest with mtry = ',
            ranked_combinations[[1,"mtry"]],
            " and ntree = ", 
            ranked_combinations[[1,"ntree"]])
    
    final_list <- list(best_model = final_model$model,
               grid_performance = ranked_combinations)
    
    return(final_list)
  }

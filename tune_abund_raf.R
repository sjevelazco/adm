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
           metrics,
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
    
    # names(hyper_combinations) <- grid[i,"comb_id"]
    # hyper_combinations <- dplyr::bind_rows(hyper_combinations, .id = "comb_id")
    
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
        mtry = ranked_combinations[[1,"mtry"]],
        ntree = ranked_combinations[[1,"ntree"]]
      )
    
    message('The best model was a Random Forest with mtry = ',
            ranked_combinations[[1,"mtry"]],
            " and ntree = ", 
            ranked_combinations[[1,"ntree"]])
    
    final_list <- append(final_model,list(ranked_combinations = ranked_combinations))
    
    return(final_list)
  }

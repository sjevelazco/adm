# source("D:/DOCUMENTOS/R/utils/para_testes_adm.R") # data, pred, response
# data = data
# response = response
# predictors = predictors
# predictors_f = NULL
# fit_formula = NULL
# partition = "kfold"
# predict_part = FALSE
# max_depth <- seq(from = 1, to = 2*length(predictors), by = 1)
# eta <- seq(from = 0.05, to = 0.2, by = 0.05)
# objective <- c("reg:squarederror","count:poisson","reg:tweedie")
# grid <- expand.grid(objective = objective, max_depth = max_depth, eta = eta)
# metrics = c("corr_spear","corr_pear")
# n_cores = 6
# verbose = FALSE
#
tune_abund_gbm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = FALSE) {
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Gradient Boosting Machines.")
      max_depth <- seq(from = 1, to = 2*length(predictors), by = 1)
      eta <- seq(from = 0.05, to = 0.2, by = 0.05)
      objective <- c("reg:squarederror","reg:logistic","count:poisson","reg:tweedie")
      nrounds <- seq(from = 250, to = 500, by = 250)
      grid <- expand.grid(objective = objective, 
                          max_depth = max_depth, 
                          eta = eta, 
                          nrounds = nrounds)
    } else {
      if (all(names(grid) %in% c("objective","max_depth", "eta","nrounds")) & length(names(grid)) == 3) {
        grid <- grid
      } else {
        stop("Grid names expected to be objective, max_depth, eta and nrounds.")
      }
    }
    
    grid$comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    hyper_combinations <- list()
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .export = c("fit_abund_gbm"), .packages = c("dplyr")) %dopar% {
      # for (i in 1:nrow(grid)) {
      # cat(i,"/",nrow(grid)) # DEBUG
      model <-
        fit_abund_gbm(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          params = list(
            max_depth = grid[i,"max_depth"],
            eta = grid[i,"eta"],
            objective = grid[i,"objective"]
          ),
          nrounds = grid[i,"nrounds"],
          verbose = 0
        )
      # hyper_combinations <- append(hyper_combinations,list(model$performance))
      # names(hyper_combinations)[i] <- grid[i,"comb_id"]
      l <- list(model$performance)
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)
    
    hyper_combinations <- dplyr::bind_rows(hyper_combinations, .id = "comb_id")
    
    ranked_combinations <- model_selection(hyper_combinations, metrics) %>%
      dplyr::left_join(grid, by = c("comb_id" = "comb_id"))
    
    for (i in 1:length(ranked_combinations)) {
      ranked_combinations[[i]] <- ranked_combinations[[i]] %>%
        dplyr::left_join(grid, by = c("comb_id" = "comb_id"))
    }
    
    # fit final model
    message("Fitting the best model...")
    # final_model <-
    #   fit_abund_raf(
    #     data = data,
    #     response = response,
    #     predictors = predictors,
    #     predictors_f = predictors_f,
    #     fit_formula = fit_formula,
    #     partition = partition,
    #     predict_part = predict_part,
    #     mtry = ranked_combinations[[1]][[1, "mtry"]],
    #     ntree = ranked_combinations[[1]][[1, "ntree"]]
    #   )
    
    message(
      "The best model was a Random Forest with mtry = ",
      ranked_combinations[[1]][[1, "mtry"]],
      " and ntree = ",
      ranked_combinations[[1]][[1, "ntree"]]
    )
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }
rm(tune_abund_gbm)


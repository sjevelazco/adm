#' Fit and validate Generalized Boosted Regression models with exploration of hyper-parameters that optimize performance
#'
#' @param data 
#' @param response 
#' @param predictors 
#' @param predictors_f 
#' @param fit_formula 
#' @param partition 
#' @param predict_part 
#' @param grid 
#' @param distribution 
#' @param metrics 
#' @param n_cores 
#' @param verbose 
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' 
#' @return
#' @export
#'
#' @examples
tune_abund_gbm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           distribution,
           metrics = NULL,
           n_cores = 1,
           verbose = FALSE) {
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    grid_dict <- list(n.trees = c(100,200,300),
                      interaction.depth = c(1,2,3),
                      n.minobsinnode = c(5,10,15),
                      shrinkage = seq(0.001,0.1,by=0.05))
    
    
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Gradient Boosting Machines.")
      grid <- expand.grid(grid_dict)
    } else if (all(names(grid) %in% names(grid_dict))) {
      user_hyper <- names(grid)[which(names(grid_dict) == names(grid))]
      default_hyper <- names(grid_dict)[which(names(grid_dict) != user_hyper)]
      
      user_list <- grid_dict[default_hyper]
      for (i in user_hyper) {
        l <- grid[[i]] %>% unique () %>%list()
        names(l) <- i
        user_list <- append(user_list,l)
      }
      
      grid <- expand.grid(user_list)
      if (all(names(grid) %in% names(grid_dict)) & length(names(grid))==4){
        message("Using provided grid.")
      }
    } else {
      stop('Grid expected to be any combination between "n.trees", "interaction.depth", "n.minobsinnode" and "shrinkage"  hyperparameters.')
    }
    
    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id,grid)
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid),.options.snow = opts, .export = c("fit_abund_gbm","adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_gbm(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          distribution = distribution,
          n.trees = grid[i,"n.trees"],
          interaction.depth = grid[i,"interaction.depth"],
          n.minobsinnode = grid[i,"n.minobsinnode"],
          shrinkage = grid[i,"shrinkage"],
          verbose = verbose
        )
      l <- list(cbind(grid[i,], model$performance))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)
    
    hyper_combinations <- lapply(hyper_combinations, function(x) dplyr::bind_rows(x)) %>% 
      dplyr::bind_rows()
      
    ranked_combinations <- model_selection(hyper_combinations, metrics)
    
    # fit final model
    message("Fitting the best model...")
    final_model <-
      fit_abund_gbm(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        distribution = distribution,
        n.trees = ranked_combinations[[1]][1,"n.trees"],
        interaction.depth = ranked_combinations[[1]][1,"interaction.depth"],
        n.minobsinnode = ranked_combinations[[1]][1,"n.minobsinnode"],
        shrinkage = ranked_combinations[[1]][1,"shrinkage"],
        verbose = verbose
      )
    
    message(
      "The best model was a GBM with n.trees = ",
      ranked_combinations[[1]][1,"n.trees"],
      ", interaction.depth = ",
      ranked_combinations[[1]][1,"interaction.depth"],
      ", n.minobsinnode = ",
      ranked_combinations[[1]][1,"n.minobsinnode"],
      " and shrinkage = ",
      ranked_combinations[[1]][1,"shrinkage"]
    )
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }



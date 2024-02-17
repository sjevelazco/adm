# data("sppabund")
# data <- sppabund
# response <- "ind_ha"
# predictors <- c("bio1", "bio12", "bio15",  "bio3",  "cfvo", "elevation","sand")
# predictors_f = NULL
# fit_formula = NULL
# partition <- ".part"
# predict_part = FALSE
# grid = NULL
# family = "poisson"
# weight_rules = c("!=0"=1)
# metrics = c("corr_spear")
# n_cores = 10
# verbose = FALSE
#
tune_abund_glm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           family,
           weight_rules = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = FALSE) {
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f)) %>%
      t() %>%
      as.vector()
    
    ziformula = c()
    for (i in 1:length(variables)) {
      pred_zero <- combn(variables,i) %>% 
        t()
      formulae <- rep(NA,nrow(pred_zero))
      for (j in 1:nrow(pred_zero)) {
        formulae[j] <- paste0("~",paste(pred_zero[j,],collapse = "+"))
      }
      ziformula <- append(ziformula,formulae)
    }
    
    ziformula <- append(ziformula,"~0")
    
    if(!(is.null(weight_rules))){
      wgts <- weight_rules
      
      resp_column <- paste0("data[,'",response,"']")
      values <- list()
      for (i in 1:length(wgts)) {
        conditions <- strsplit(names(wgts)[i],"\\s*&\\s*")[[1]]
        if(length(conditions)==2){
          logical_test <- paste0(paste0(resp_column,conditions), collapse = "&")
        } else if (length(conditions)==1){
          logical_test <- paste0(resp_column,names(wgts)[i])
        }
        result <- ifelse(eval(parse(text = logical_test)),wgts[i],0) %>% 
          as.vector()
        values <- append(values,list(result))
      }
      
      weights_obs <- Reduce("+",values)
      weights_obs <- replace(weights_obs,weights_obs==0,1)
      weights_obs <- weights_obs/sum(weights_obs)
      
      weighted <- c(TRUE,FALSE)
    } else {
      weights_obs <- NULL
      weighted <- c(FALSE)
    }
    
    
    grid_dict <- list(ziformula = ziformula,
                      weighted = weighted)
    
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for GLM.")
      
      grid <- expand.grid(grid_dict)
    } else if (all(names(grid) %in% c("ziformula","weighted"))) {
      user_hyper <- names(grid)[which(names(grid_dict) == names(grid))]
      default_hyper <- names(grid_dict)[which(names(grid_dict) != user_hyper)]
      
      user_list <- grid_dict[default_hyper]
      for (i in user_hyper) {
        l <- grid[[i]] %>% unique () %>%list()
        names(l) <- i
        user_list <- append(user_list,l)
      }
      
      grid <- expand.grid(user_list)
      if (all(names(grid) %in% c("weighted","ziformula")) & length(names(grid))==2){
        message("Using provided grid.")
      } else {
        stop('Grid expected to be any combination between "ziformula" and "weighted" hyperparameters.')
      }
    } else {
      stop('Grid expected to be any combination between "ziformula" and "weighted" hyperparameters.')
    }
    
    
    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id,grid)
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .export = c("fit_abund_glm","adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_glm(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          family = family,
          weights_obs = if(grid[i,"weighted"]==FALSE){NULL}else{weights_obs},
          ziformula = formula(as.character(grid[i,"ziformula"]))
        )
      l <- list(cbind(grid[i,], model$performance))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)
    
    hyper_combinations <- lapply(hyper_combinations, function(x) bind_rows(x)) %>% 
      bind_rows()
    
    ranked_combinations <- model_selection(hyper_combinations, metrics)
    
    # fit final model
    message("Fitting the best model...")
    final_model <-
      fit_abund_glm(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        family = family,
        weights_obs = if(ranked_combinations[[1]][1,"weighted"]==FALSE){NULL}else{weights_obs},
        ziformula = formula(as.character(ranked_combinations[[1]][1,"ziformula"]))
      )

    message(
      "The best model was a GLM with:", 
      "\n ziformula = ",
      ranked_combinations[[1]][1,"ziformula"],
      "\n weighted = ",
      ranked_combinations[[1]][1,"weighted"]
    )
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }

#' Fit and validate Deep Neural Network model with exploration of hyper-parameters that optimize performance
#' 
#' @param data 
#' @param response 
#' @param predictors 
#' @param predictors_f 
#' @param fit_formula 
#' @param partition 
#' @param predict_part 
#' @param grid 
#' @param architectures 
#' @param metrics 
#' @param n_cores 
#' @param verbose 
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stringr str_extract_all
#' 
#' @return
#' @export
#'
#' @examples
tune_abund_dnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           architectures = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = FALSE) {
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    # architectures
    if (all(architectures != "fit_intern")){
      if (is.null(architectures)){
        message("Architectures not provided. Using the default set for Deep Neural Networks.")
        architectures <- generate_arch_list(number_of_features = length(predictors)+length(predictors_f),
                                            number_of_outputs = length(response))
        arch_list <- architectures$arch_list
        arch_dict <- architectures$arch_dict
      } else if (!all(names(architectures)%in%c("arch_list","arch_dict"))){
        stop("architectures expected to be a list with two other lists, arch_list and arch_dict, or 'fit_intern'.")
      } else {
        arch_list <- architectures$arch_list
        arch_dict <- architectures$arch_dict
        if (!all(sapply(arch_list, class) == c("neural_net","nn_module","nn_module_generator"))){
          stop('Expected a "neural_net","nn_module","nn_module_generator" objects in arch_list.')
        } else {
          message('Using provided architectures.')
        }
      }
    } else {
      message("Using the fit_abund_dnn() intern DNN architecture.")
      arch_list <- list("fit_intern"=NULL)
      arch_dict <- list("fit_intern"=NULL)
    }

    archs <- names(arch_list)
    
    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Deep Neural Networks.")
      batch_size <- 2^seq(4,6)
      #n_epochs <- seq(10,20,by=10)
      n_epochs <- 10
      learning_rate <- seq(from = 0.1, to = 0.2, by = 0.1)
      grid <- expand.grid(arch = archs,
                          batch_size = batch_size, 
                          n_epochs = n_epochs, 
                          learning_rate = learning_rate)
    } else {
      if (all(names(grid) %in% c("batch_size","n_epochs", "learning_rate")) & length(names(grid)) == 3) {
        batch_size <- unique(grid[,"batch_size"])
        n_epochs <- unique(grid[,"n_epochs"])
        learning_rate <- unique(grid[,"learning_rate"])
        grid <- expand.grid(arch = archs,
                            batch_size = batch_size, 
                            n_epochs = n_epochs, 
                            learning_rate = learning_rate)
      } else {
        stop("Grid names expected to be batch_size, n_epochs and learning_rate")
      }
    }
    
    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id,grid)
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .export = c("fit_abund_dnn","adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_dnn(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          learning_rate = grid[i,"learning_rate"],
          n_epochs = grid[i,"n_epochs"],
          batch_size = grid[i,"batch_size"],
          custom_architecture = arch_list[[grid[i,"arch"]]]
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
      fit_abund_dnn(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        learning_rate = ranked_combinations[[1]][1,"learning_rate"],
        n_epochs = ranked_combinations[[1]][1,"n_epochs"],
        batch_size = ranked_combinations[[1]][1,"batch_size"],
        custom_architecture = arch_list[[ranked_combinations[[1]][1,"arch"]]]
      )
    
    arch_indexes <- stringr::str_extract_all(ranked_combinations[[1]][1,"arch"], "\\d+")

    message(
      "The best model was a DNN with learning_rate = ",
      ranked_combinations[[1]][1,"learning_rate"],
      ", n_epochs = ",
      ranked_combinations[[1]][1,"n_epochs"],
      ", batch_size = ",
      ranked_combinations[[1]][1,"batch_size"]
    )
    
    if (ranked_combinations[[1]][1,"arch"] != "fit_intern"){
      n_layers <- as.numeric(arch_indexes[[1]][1])
      n_comb <- as.numeric(arch_indexes[[1]][2])
      
      structure <- arch_dict[[n_layers]][,n_comb]
      message("Used a ", n_layers, " hidden layers DNN structured as ",
              paste0(structure, collapse = "-"))
    }
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }



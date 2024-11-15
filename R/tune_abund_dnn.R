#' Fit and validate Deep Neural Network model with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with "batch_size", "n_epochs", "learning_rate" as columns and its values combinations as rows. If no grid is provided, function will create a default grid combining the next hyperparameters:
#' batch_size = 2^seq(4, 6)
#' n_epochs = 10
#' learning_rate = seq(from = 0.1, to = 0.2, by = 0.1)
#' validation_patience = 2
#' fitting_patience = 5. In case one or more hyperparameters are provided, the function will complete the grid with the default values.
#' @param architectures list or character. A list object containing a list of architectures (nn_modules_generators from torch), called "arch_list", and a list of matrices describing each architecture, called ("arch_dict"); use generate_arch_list function to create it. It's also possible to use "fit_intern", what will construct the default neural network architecture of fit_abund_dnn. If NULL, a list of architectures will be generated. Default NULL
#' @param metrics character. Vector with one or more metrics from c("corr_spear","corr_pear","mae","pdisp","inter","slope").
#' @param n_cores numeric. Number of cores used in parallel processing.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stringr str_extract_all
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "luz_module_fitted" object from luz (torch framework). This object can be used to predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: A tibble with selected model's performance metrics calculated in adm_eval.
#' \item performance_part: A tibble with performance metrics for each test partition.
#' \item predicted_part: A tibble with predicted abundance for each test partition.
#' \item optimal_combination: A tibble with the selected hyperparameter combination and its performance.
#' \item all_combinations: A tibble with all hyperparameters combinations and its performance.
#' \item selected_arch: A numeric vector describing the selected architecture layers.
#' }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' 
#' # Database with species abundance and x and y coordinates
#' data("sppabund")
#' 
#' # Select data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(-.part2, -.part3)
#'   
#' # Explore response variables
#' some_sp$ind_ha %>% range()
#' some_sp$ind_ha %>% hist()
#' 
#' # Here we balance number of absences
#' some_sp <-
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
#' 
#' # Generate some architectures
#' many_archs <- generate_arch_list(
#'   type = "dnn",
#'   number_of_features = 3,
#'   number_of_outputs = 1,
#'   n_layers = c(2,3),
#'   n_neurons = c(6, 12, 16)
#' ) %>% select_arch_list(
#'   type = c("dnn"),
#'   method = "percentile",
#'   n_samples = 1,
#'   min_max = TRUE
#' )
#' 
#' # Create a grid
#' # Obs.: the grid is tested with every architecture, thus it can get very large.
#' dnn_grid <- expand.grid(
#'   learning_rate = c(0.01, 0.005),
#'   n_epochs = c(50,100),
#'   batch_size = c(32),
#'   validation_patience = c(2,4),
#'   fitting_patience = c(2,4)
#' )
#' 
#' # Tune a DNN model
#' tuned_dnn <- tune_abund_dnn(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = dnn_grid,
#'   architectures = many_archs,
#'   n_cores = 3
#' )
#' 
#' tuned_dnn
#' 
#' # It is also possible to use a only one architecture 
#' one_arch <- generate_dnn_architecture(
#'   number_of_features = 3,
#'   number_of_outputs = 1,
#'   number_of_hidden_layers = 3,
#'   hidden_layers_size = c(8, 16, 8),
#'   batch_norm = TRUE
#' )
#' 
#' tuned_dnn_2 <- tune_abund_dnn(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = dnn_grid,
#'   architectures = one_arch,
#'   n_cores = 3
#' )
#' 
#' tuned_dnn_2
#' }
tune_abund_dnn <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           architectures = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = TRUE) {
    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }

    # architectures
    if (is.list(architectures)){
      # check if it is from generate_arch_list or generate_dnn_architecture
      if (all(names(architectures) %in% c("net","arch","arch_dict"))){
        # generated with generate_dnn_architecture
        arch_dict <- architectures$arch_dict
        arch_list <- list(architectures$net)
        names(arch_list) <- paste0("arch-",arch_dict %>% names() %>% stringr::str_replace_all("[^0-9]", ""),"-1")
      } else if (all(names(architectures) %in% c("arch_list","arch_dict","changes"))) {
        # generated with generate_arch_list
        arch_list <- architectures$arch_list
        arch_dict <- architectures$arch_dict
      }
    } else if (is.null(architectures) | all(architectures == "fit_intern")) {
      # use the default Deep Neural Network architecture
      message("Using the fit_abund_dnn() intern DNN architecture.")
      arch_list <- list("fit_intern" = NULL)
      arch_dict <- list("fit_intern" = NULL)
    }
    
    if (all(names(arch_list) != c("fit_intern"))){
      if (!all(sapply(arch_list, class) == c("neural_net", "nn_module", "nn_module_generator"))) {
        stop('Expected "neural_net", "nn_module", "nn_module_generator" objects in arch_list.
      Please, use generate_arch_list or generate_dnn_architecture outputs.')
      } else {
        message("Using provided architectures.")
      }
    }

    # making grid
    grid_dict <- list(
      learning_rate = c(0.01, 0.005),
      n_epochs = c(100,200),
      batch_size = c(16,32),
      validation_patience = c(2,4),
      fitting_patience = c(2,4)
    )
    
    archs <- names(arch_list)
    
    # Check hyperparameters names
    nms_grid <- names(grid)
    nms_hypers <- names(grid_dict)
    
    if (!all(nms_grid %in% nms_hypers)) {
      stop(
        paste(paste(nms_grid[!nms_grid %in% nms_hypers], collapse = ", "), " is not hyperparameters\n"),
        "Grid expected to be any combination between ", paste(nms_hypers, collapse = ", ")
      )
    }
    
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Support Vector Machines.")
      
      grid_dict <- append(grid_dict, list(arch = archs))
      
      grid <- expand.grid(grid_dict)
    } else if (all(nms_hypers %in% nms_grid)) {
      user_list <- list()
      for (i in nms_grid) {
        l <- grid[[i]] %>%
          unique() %>%
          list()
        names(l) <- i
        user_list <- append(user_list, l)
      }
      
      user_list <- append(user_list, list(arch = archs))
      
      grid <- expand.grid(user_list)
      
      message("Using provided grid.")
    } else if (any(!nms_hypers %in% nms_grid)) {
      message(
        "Adding default hyperparameter for: ",
        paste(names(grid_dict)[!names(grid_dict) %in% nms_grid], collapse = ", ")
      )
      
      user_hyper <- names(grid)[which(names(grid) %in% names(grid_dict))]
      default_hyper <- names(grid_dict)[which(!names(grid_dict) %in% user_hyper)]
      
      user_list <- grid_dict[default_hyper]
      for (i in user_hyper) {
        l <- grid[[i]] %>%
          unique() %>%
          list()
        names(l) <- i
        user_list <- append(user_list, l)
      }
      
      user_list <- append(user_list, list(arch = archs))
      
      grid <- expand.grid(user_list)
    }

    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)
    message(paste0("Testing ",nrow(grid)," combinations."))

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_dnn", "adm_eval", "adapt_df"), .packages = c("dplyr", "torch")) %dopar% {
      model <-
        fit_abund_dnn(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          partition = partition,
          predict_part = predict_part,
          learning_rate = grid[i, "learning_rate"],
          n_epochs = grid[i, "n_epochs"],
          batch_size = grid[i, "batch_size"],
          custom_architecture = arch_list[[grid[i, "arch"]]],
          validation_patience = grid[i, "validation_patience"],
          fitting_patience = grid[i, "fitting_patience"], 
          verbose = verbose
        )
      l <- list(cbind(grid[i, ], model$performance))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)

    hyper_combinations <- lapply(hyper_combinations, function(x) dplyr::bind_rows(x)) %>%
      dplyr::bind_rows()

    for (i in 1:ncol(hyper_combinations)) {
      if (all(is.na(hyper_combinations[[i]]))) {
        stop(paste0("The net was unable to fit the data. Try changing the hyperparameters."))
      }
    }

    ranked_combinations <- model_selection(hyper_combinations, metrics)

    # fit final model
    message("\nFitting the best model...")
    final_model <-
      fit_abund_dnn(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        partition = partition,
        predict_part = predict_part,
        learning_rate = ranked_combinations[[1]][1, "learning_rate"],
        n_epochs = ranked_combinations[[1]][1, "n_epochs"],
        batch_size = ranked_combinations[[1]][1, "batch_size"],
        custom_architecture = arch_list[[ranked_combinations[[1]][1, "arch"]]],
        validation_patience = ranked_combinations[[1]][1, "validation_patience"],
        fitting_patience = ranked_combinations[[1]][1, "fitting_patience"],
        verbose = verbose
      )

    arch_indexes <- stringr::str_extract_all(ranked_combinations[[1]][1, "arch"], "\\d+")
    n_layers <- arch_indexes[[1]][[1]]
    layer_index <- paste0(n_layers, "_layer_net")
    size_index <- arch_indexes[[1]][[2]]

    message(
      "The best model was achieved with: \n learning_rate = ",
      ranked_combinations[[1]][1, "learning_rate"],
      "\n n_epochs = ",
      ranked_combinations[[1]][1, "n_epochs"],
      "\n patience = ",
      ranked_combinations[[1]][1, "validation_patience"],
      " and ",
      ranked_combinations[[1]][1, "fitting_patience"],
      "\n batch_size = ",
      ranked_combinations[[1]][1, "batch_size"],
      "\n arch = ",
      n_layers, " layers with ",
      paste(arch_dict[[layer_index]][, size_index %>% as.numeric()], collapse = "->"), " neurons "
    )

    selected_arch <- ranked_combinations[[1]][1, "arch"] %>%
      as.character()
    selected_arch <- gsub("arch-", "", selected_arch)
    selected_arch <- paste0(substr(selected_arch, 1, nchar(selected_arch) - 2), "_layer_net")
    n_comb <- as.numeric(arch_indexes[[1]][2])

    final_list <- c(final_model, ranked_combinations, list("selected_arch" = arch_dict[[selected_arch]][, n_comb]))

    return(final_list)
  }

#' Fit and validate Deep Neural Network model with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with "batch_size", "n_epochs", "learning_rate" as columns and its values combinations as rows.
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
           verbose = FALSE) {
    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }

    # architectures
    if (all(architectures != "fit_intern")) {
      if (is.null(architectures)) {
        message("Architectures not provided. Using the default set for Deep Neural Networks.")
        architectures <- generate_arch_list(
          number_of_features = length(predictors) + length(predictors_f),
          number_of_outputs = length(response)
        )
        arch_list <- architectures$arch_list
        arch_dict <- architectures$arch_dict
      } else if (!all(names(architectures) %in% c("arch_list", "arch_dict", "changes"))) {
        stop("architectures expected to be a list with two other lists, arch_list and arch_dict, or 'fit_intern'.")
      } else {
        arch_list <- architectures$arch_list
        arch_dict <- architectures$arch_dict
        if (!all(sapply(arch_list, class) == c("neural_net", "nn_module", "nn_module_generator"))) {
          stop('Expected "neural_net", "nn_module", "nn_module_generator" objects in arch_list.')
        } else {
          message("Using provided architectures.")
        }
      }
    } else {
      message("Using the fit_abund_dnn() intern DNN architecture.")
      arch_list <- list("fit_intern" = NULL)
      arch_dict <- list("fit_intern" = NULL)
    }

    archs <- names(arch_list)

    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Deep Neural Networks.")
      batch_size <- 2^seq(4, 6)
      # n_epochs <- seq(10,20,by=10)
      n_epochs <- 10
      learning_rate <- seq(from = 0.1, to = 0.2, by = 0.1)
      validation_patience <- 2
      fitting_patience <- 5
      grid <- expand.grid(
        arch = archs,
        batch_size = batch_size,
        n_epochs = n_epochs,
        learning_rate = learning_rate,
        validation_patience = validation_patience,
        fitting_patience = fitting_patience
      )
    } else {
      if (all(names(grid) %in% c("batch_size", "n_epochs", "learning_rate", "validation_patience","fitting_patience")) & length(names(grid)) == 5) {
        batch_size <- unique(grid[, "batch_size"])
        n_epochs <- unique(grid[, "n_epochs"])
        learning_rate <- unique(grid[, "learning_rate"])
        validation_patience <- unique(grid[, "validation_patience"])
        fitting_patience <- unique(grid[, "fitting_patience"])
        grid <- expand.grid(
          arch = archs,
          batch_size = batch_size,
          n_epochs = n_epochs,
          learning_rate = learning_rate,
          validation_patience = validation_patience,
          fitting_patience = fitting_patience
        )
      } else {
        stop("Grid names expected to be batch_size, n_epochs, learning_rate, validation_patience and fitting_patience")
      }
    }

    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_dnn", "adm_eval", "nnf_dropout","adapt_df"), .packages = c("dplyr", "torch")) %dopar% {
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
          fitting_patience = grid[i, "fitting_patience"]
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
    message("Fitting the best model...")
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
        fitting_patience = ranked_combinations[[1]][1, "fitting_patience"]
      )

    arch_indexes <- stringr::str_extract_all(ranked_combinations[[1]][1, "arch"], "\\d+")
    n_layers <- arch_indexes[[1]][[1]]
    layer_index <- paste0(n_layers, "_layer_net")
    size_index <- arch_indexes[[1]][[2]]

    message(
      "The best model was a DNN with: \n learning_rate = ",
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

    # if (ranked_combinations[[1]][1, "arch"] != "fit_intern") {
    #   n_layers <- as.numeric(arch_indexes[[1]][1])
    #   n_comb <- as.numeric(arch_indexes[[1]][2])
    #
    #   structure <- arch_dict[[n_layers]][, n_comb]
    #   message(
    #     "Used a ", n_layers, " hidden layers DNN structured as ",
    #     paste0(structure, collapse = "-")
    #   )
    # }

    selected_arch <- ranked_combinations[[1]][1, "arch"] %>%
      as.character()
    selected_arch <- gsub("arch-", "", selected_arch)
    selected_arch <- paste0(substr(selected_arch, 1, nchar(selected_arch) - 2), "_layer_net")
    n_comb <- as.numeric(arch_indexes[[1]][2])

    final_list <- c(final_model, ranked_combinations, list("selected_arch" = arch_dict[[selected_arch]][, n_comb]))

    return(final_list)
  }

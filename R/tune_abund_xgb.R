#' Fit and validate Extreme Gradient Boosting models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with "n.trees", "interaction.depth", "n.minobsinnode" and "shrinkage" as columns and its values combinations as rows. If no grid is provided, function will
#' create a default grid combining the next hyperparameters:
#' nrounds = c(100, 200, 300),
#' max_depth = c(4, 6, 8),
#' eta = c(0.2, 0.4, 0.5),
#' gamma = c(1, 5, 10),
#' colsample_bytree = c(0.5, 1, 2),
#' min_child_weight = c(0.5, 1, 2),
#' subsample = c(0.5, 0.75, 1).
#' In case one or more hyperparameters are provided, the function will complete
#' the grid with the default values.
#' @param objective character. The learning task and the corresponding learning objective. Default is "reg:squarederror", regression with squared loss.
#' @param metrics character. Vector with one or more metrics from c("corr_spear","corr_pear","mae","pdisp","inter","slope").
#' @param n_cores numeric. Number of cores used in parallel processing.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows as_tibble
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#'
#' @return
#' A list object with:
#' \itemize{
#' \item model: A "xgb.Booster" object from xgboost package. This object can be used to predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: A tibble with selected model's performance metrics calculated in adm_eval.
#' \item performance_part: A tibble with performance metrics for each test partition.
#' \item predicted_part: A tibble with predicted abundance for each test partition.
#' \item optimal_combination: A tibble with the selected hyperparameter combination and its performance.
#' \item all_combinations: A tibble with all hyperparameters combinations and its performance.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' # Database with species abundance and x and y coordinates
#' data("sppabund")
#' # Select data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species two") %>%
#'   dplyr::select(-.part2, -.part3)
#' # Explore response variables
#' some_sp$ind_ha %>% range()
#' some_sp$ind_ha %>% hist()
#' # Here we balance number of absences
#' some_sp <-
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
#' # Create a grid
#' xgb_grid <- expand.grid(
#'   nrounds = c(100, 300),
#'   max_depth = c(4, 6, 8),
#'   eta = c(0.2, 0.5),
#'   gamma = c(1, 5, 10),
#'   colsample_bytree = c(0.5, 1),
#'   min_child_weight = c(0.5, 1, 2),
#'   subsample = c(0.5, 1),
#'   stringsAsFactors = FALSE
#' )
#' # Tune a XGB model
#' tuned_xgb <- tune_abund_xgb(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = xgb_grid,
#'   objective = "reg:squarederror",
#'   n_cores = 3
#' )
#' tuned_xgb
#' }
tune_abund_xgb <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           objective = "reg:squarederror",
           metrics = NULL,
           n_cores = 1,
           verbose = TRUE) {
    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }

    grid_dict <- list(
      nrounds = c(100, 200, 300),
      max_depth = c(4, 6, 8),
      eta = c(0.2, 0.4, 0.5),
      gamma = c(1, 5, 10),
      colsample_bytree = c(0.5, 1, 2),
      min_child_weight = c(0.5, 1, 2),
      subsample = c(0.5, 0.75, 1)
    )


    # making grid
    grid_dict <- list(
      nrounds = c(100, 300),
      max_depth = c(4, 6, 8),
      eta = c(0.2, 0.5),
      gamma = c(1, 5, 10),
      colsample_bytree = c(0.5, 1),
      min_child_weight = c(0.5, 1, 2),
      subsample = c(0.5, 1)
    )

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
      message("Grid not provided. Using the default one for Extreme Gradient Boosting.")
      grid <- expand.grid(grid_dict, stringsAsFactors = FALSE)
    } else if (all(nms_hypers %in% nms_grid)) {
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

      grid <- expand.grid(user_list, stringsAsFactors = FALSE)
    }

    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    # doSNOW::registerDoSNOW(cl)
    # pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    # progress <- function(n) utils::setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .export = c("fit_abund_xgb", "adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_xgb(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          partition = partition,
          predict_part = predict_part,
          max_depth = grid[i, "max_depth"],
          eta = grid[i, "eta"],
          gamma = grid[i, "gamma"],
          colsample_bytree = grid[i, "colsample_bytree"],
          min_child_weight = grid[i, "min_child_weight"],
          subsample = grid[i, "subsample"],
          objective = objective,
          nrounds = grid[i, "nrounds"],
          verbose = verbose
        )
      l <- list(cbind(grid[i, ], model$performance))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)

    hyper_combinations <- lapply(hyper_combinations, function(x) dplyr::bind_rows(x)) %>%
      dplyr::bind_rows()

    ranked_combinations <- model_selection(hyper_combinations, metrics)

    # fit final model
    message("\nFitting the best model...")
    final_model <-
      fit_abund_xgb(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        partition = partition,
        predict_part = predict_part,
        max_depth = ranked_combinations[[1]][1, "max_depth"],
        eta = ranked_combinations[[1]][1, "eta"],
        gamma = ranked_combinations[[1]][1, "gamma"],
        colsample_bytree = ranked_combinations[[1]][1, "colsample_bytree"],
        min_child_weight = ranked_combinations[[1]][1, "min_child_weight"],
        subsample = ranked_combinations[[1]][1, "subsample"],
        objective = objective,
        nrounds = ranked_combinations[[1]][1, "nrounds"],
        verbose = verbose
      )

    message(
      "The best model was achieved with: \n max_depth = ",
      ranked_combinations[[1]][1, "max_depth"],
      "\n eta = ",
      ranked_combinations[[1]][1, "eta"],
      "\n gamma = ",
      ranked_combinations[[1]][1, "gamma"],
      "\n colsample_bytree = ",
      ranked_combinations[[1]][1, "colsample_bytree"],
      "\n min_child_weight = ",
      ranked_combinations[[1]][1, "min_child_weight"],
      "\n subsample = ",
      ranked_combinations[[1]][1, "subsample"],
      "\n nrounds = ",
      ranked_combinations[[1]][1, "nrounds"]
    )

    final_list <- c(final_model, ranked_combinations)

    # Standardize output list
    for (i in 2:length(final_list)) {
      if (!class(final_list[[i]])[1] == "tbl_df") {
        final_list[[i]] <- dplyr::as_tibble(final_list[[i]])
      }
    }

    return(final_list)
  }

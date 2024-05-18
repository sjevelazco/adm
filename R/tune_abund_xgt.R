#' Fit and validate Extreme Gradient Boosting models with exploration of hyper-parameters that optimize performance
#'
#' @param data
#' @param response
#' @param predictors
#' @param predictors_f
#' @param fit_formula
#' @param partition
#' @param predict_part
#' @param grid
#' @param objective
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
tune_abund_xgt <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           objective,
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
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for eXtreme Gradient Boosting.")
      grid <- expand.grid(grid_dict)
    } else if (all(names(grid) %in% names(grid_dict))) {
      user_hyper <- names(grid)[which(names(grid_dict) == names(grid))]
      default_hyper <- names(grid_dict)[which(names(grid_dict) != user_hyper)]

      user_list <- grid_dict[default_hyper]
      for (i in user_hyper) {
        l <- grid[[i]] %>%
          unique() %>%
          list()
        names(l) <- i
        user_list <- append(user_list, l)
      }

      grid <- expand.grid(user_list)
      if (all(names(grid) %in% names(grid_dict)) & length(names(grid)) == 7) {
        message("Using provided grid.")
      }
    } else {
      stop('Grid expected to be any combination between "n.trees", "interaction.depth", "n.minobsinnode" and "shrinkage"  hyperparameters.')
    }

    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_xgt", "adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_xgt(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
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
    message("Fitting the best model...")
    final_model <-
      fit_abund_xgt(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
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
      "The best model was a XGBoost with: \n max_depth = ",
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

    return(final_list)
  }

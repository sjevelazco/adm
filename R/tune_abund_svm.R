#' Fit and validate Support Vector Machine models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param grid tibble or data.frame. A dataframe with "kernel", "sigma", "C" as columns and its values combinations as rows.
#' @param metrics character. Vector with one or more metrics from c("corr_spear","corr_pear","mae","pdisp","inter","slope").
#' @param n_cores numeric. Number of cores used in parallel processing.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom dplyr bind_rows
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "ksvm" class object from kernlab package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: A tibble with selected model's performance metrics calculated in adm_eval.
#' \item performance_part: A tibble with performance metrics for each test partition.
#' \item predicted_part: A tibble with predicted abundance for each test partition.
#' \item optimal_combination: A tibble with the selected hyperparameter combination and its performance.
#' \item all_combinations: A tibble with all hyperparameters combinations and its performance.
#' }
#' @export
#'
#' @examples
tune_abund_svm <-
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

    grid_dict <- list(
      C = seq(0.2, 1, by = 0.2),
      sigma = "automatic",
      kernel = c("rbfdot", "laplacedot")
    )

    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Support Vector Machines.")
      grid <- expand.grid(grid_dict)
    } else if (all(names(grid) %in% c("C", "sigma", "kernel"))) {
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
      if (all(names(grid) %in% c("C", "sigma", "kernel")) & length(names(grid)) == 3) {
        message("Using provided grid.")
      }
    } else {
      stop('Grid expected to be any combination between "c", "sigma" and "kernel" hyperparameters.')
    }


    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)
    grid[["kernel"]] <- as.character(grid[["kernel"]])
    grid[["sigma"]] <- as.character(grid[["sigma"]])

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_svm", "adm_eval"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_svm(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          sigma = grid[[i, "sigma"]],
          kernel = grid[[i, "kernel"]],
          C = grid[[i, "C"]]
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
      fit_abund_svm(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        sigma = ranked_combinations[[1]][[1, "sigma"]],
        kernel = ranked_combinations[[1]][[1, "kernel"]],
        C = ranked_combinations[[1]][[1, "C"]]
      )

    message(
      "The best model was a Support Vector Machine with sigma = ",
      ranked_combinations[[1]][[1, "sigma"]],
      ", kernel = ",
      ranked_combinations[[1]][[1, "kernel"]],
      " and C = ",
      ranked_combinations[[1]][[1, "C"]]
    )

    final_list <- c(final_model, ranked_combinations)

    return(final_list)
  }

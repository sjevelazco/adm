#' Fit and validate Support Vector Machine models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param grid tibble or data.frame. A dataframe with "kernel", "sigma", "C" as columns and
#' its values combinations as rows. If now grid is provided, funcion will create a default grid combining
#' the next hyperparameters:
#' C = seq(0.2, 1, by = 0.2),
#' sigma = "automatic",
#' kernel = c("rbfdot", "laplacedot").
#' In case one or more hyperparameters are provided, the function will complete the grid with the default values.
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
#' # Create a grid
#' svm_grid <- expand.grid(
#'   sigma = "automatic",
#'   C = c(0.5, 1, 2),
#'   kernel = c("rbfdot", "laplacedot"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Tune a SVM model
#' tuned_svm <- tune_abund_svm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = svm_grid,
#'   n_cores = 3
#' )
#'
#' tuned_svm
#' }
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
           verbose = TRUE) {
    # Check metrics
    check_metrics(metrics)

    # making grid
    grid_dict <- list(
      C = seq(0.2, 1, by = 0.2),
      sigma = "automatic",
      kernel = c("rbfdot", "laplacedot")
    )

    # Check hyperparameters names
    grid <- build_search_grid(grid, grid_dict)

    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id, grid)
    grid[["kernel"]] <- as.character(grid[["kernel"]])
    grid[["sigma"]] <- as.character(grid[["sigma"]])

    # looping the grid
    message("Searching for optimal hyperparameters...")

    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    # doSNOW::registerDoSNOW(cl)
    # pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    # progress <- function(n) utils::setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)

    on.exit(
      {
        tryCatch(
          {
            parallel::stopCluster(cl)
          },
          error = function(e) {}
        )
      },
      add = T
    )
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .export = c("fit_abund_svm", "adm_eval"), .packages = c("dplyr")) %dopar% {
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
          C = grid[[i, "C"]],
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
        C = ranked_combinations[[1]][[1, "C"]],
        verbose = verbose
      )

    message(
      "The best model was achieved with: \n sigma = ",
      ranked_combinations[[1]][[1, "sigma"]],
      ", kernel = ",
      ranked_combinations[[1]][[1, "kernel"]],
      " and C = ",
      ranked_combinations[[1]][[1, "C"]]
    )

    final_list <- c(final_model, ranked_combinations)

    # # Standardize output list
    # for (i in 2:length(final_list)) {
    #   if (!class(final_list[[i]])[1] == "tbl_df") {
    #     final_list[[i]] <- dplyr::as_tibble(final_list[[i]])
    #   }
    # }

    return(final_list)
  }

#' Fit and validate Artificial Neural Networks models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with "mtry" and "ntree" as columns and its values combinations as rows.
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
#' \item model: A "nnet" class object from nnet package. This object can be used for predicting.
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
#'
#' data("sppabund")
#'
#' some_sp <- sppabund %>%
#'   filter(species == "Species two")
#'
#' net_grid <- expand.grid(
#'   size = seq(from = 8, to = 32, by = 6),
#'   decay = seq(from = 0, to = 0.4, by = 0.01)
#' )
#'
#' tuned_net <- tune_abund_net(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = net_grid,
#'   n_cores = 3
#' )
#'
#' tuned_net$model
#' tuned_net$performance
#' tuned_net$optimal_combination
#' tuned_net$all_combinations
#' }
tune_abund_net <-
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
    i <- NULL

    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }

    # making grid
    if (is.null(grid)) {
      message("Grid not provided. Using the default one for Shallow Neural Networks.")
      size <- seq(from = length(c(predictors, predictors_f)), to = 50, by = 2)
      decay <- seq(from = 0, to = 0.9, by = 0.1)
      grid <- expand.grid(size = size, decay = decay)
    } else {
      if (all(names(grid) %in% c("size", "decay")) & length(names(grid)) == 2) {
        grid <- grid
      } else {
        stop("Grid names expected to be size and decay.")
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

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_net", "adm_eval", "adapt_df"), .packages = c("dplyr")) %dopar% {
      model <-
        fit_abund_net(
          data = data,
          response = response,
          predictors = predictors,
          predictors_f = predictors_f,
          fit_formula = fit_formula,
          partition = partition,
          predict_part = predict_part,
          size = grid[i, "size"],
          decay = grid[i, "decay"]
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
      fit_abund_net(
        data = data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        size = ranked_combinations[[1]][[1, "size"]],
        decay = ranked_combinations[[1]][[1, "decay"]]
      )

    message(
      "The best model was a Artificial Neural Network with size = ",
      ranked_combinations[[1]][[1, "size"]],
      " and decay = ",
      ranked_combinations[[1]][[1, "decay"]]
    )

    final_list <- c(final_model, ranked_combinations)

    return(final_list)
  }

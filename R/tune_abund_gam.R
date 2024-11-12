#' Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param sigma_formula formula. formula for fitting a model to the nu parameter. Usage sigma_formula = ~ precipt + temp
#' @param nu_formula formula. formula for fitting a model to the nu parameter. Usage nu_formula = ~ precipt + temp
#' @param tau_formula formula. formula for fitting a model to the tau parameter. Usage tau_formula = ~ precipt + temp
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with 'distribution', 'inter' as columns and its values combinations as rows. If no grid is provided, function will create a default grid combining the next hyperparameters: distribution = all families selected by family_selector, inter = "automatic". In case one or more hyperparameters are provided, the function will complete the grid with the default values.
#'
#' @param metrics character. Vector with one or more metrics from c("corr_spear","corr_pear","mae","pdisp","inter","slope")
#' @param n_cores numeric. Number of cores used in parallel processing.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom dplyr bind_rows left_join select
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats na.omit
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gamlss" object from gamlss package. This object can be used to predicting.
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
#' require(gamlss)
#' 
#' # Database with species abundance and x and y coordinates
#' data("sppabund")
#' # Select data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(-.part2, -.part3)
#' # Explore response variables
#' some_sp$ind_ha %>% range()
#' some_sp$ind_ha %>% hist()
#' # Here we balance number of absences
#' some_sp <-
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
#' # Explore different family distributions
#' suitable_distributions <- family_selector(data = some_sp, response = "ind_ha")
#' suitable_distributions
#' # Create a grid
#' gam_grid <- expand.grid(
#'   inter = "automatic",
#'   distribution = suitable_distributions$family_call
#' )
#' # Tune a GAM model
#' tuned_gam <- tune_abund_gam(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   fit_formula = formula("ind_ha ~ bio12 + elevation + sand + eco"),
#'   sigma_formula = formula("ind_ha ~ bio12 + elevation"),
#'   nu_formula = formula("ind_ha ~ bio12 + elevation"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   predict_part = TRUE,
#'   metrics = c("corr_pear", "mae"),
#'   grid = gam_grid,
#'   n_cores = 3
#' )
#' 
#' tuned_gam
#' }
tune_abund_gam <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           sigma_formula = ~1,
           nu_formula = ~1,
           tau_formula = ~1,
           partition,
           predict_part = FALSE,
           grid = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = TRUE) {
    . <- discrete <- i <- performance <- NULL

    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    # making grid
    suitable_distributions <- family_selector(data, response)
    
    grid_dict <- list(
      inter = "automatic",
      distribution = suitable_distributions$family_call
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
      message("Grid not provided. Using the default one for Generalized Additive Models.")
      grid <- expand.grid(grid_dict)
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
      
      grid <- expand.grid(user_list)
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

    hyper_combinations <- foreach::foreach(
      i = 1:nrow(grid),
      .options.snow = opts,
      .export = c("fit_abund_gam", "adm_eval"),
      .packages = c("dplyr", "gamlss")
    ) %dopar% {
      data_fam <- data
      # if (grid[i, "discrete"] == 1) {
      #   data_fam[, response] <- round(data[, response])
      # }

      model <- tryCatch(
        {
          model <-
            fit_abund_gam(
              data = data_fam,
              response = response,
              predictors = predictors,
              predictors_f = predictors_f,
              fit_formula = fit_formula,
              partition = partition,
              predict_part = predict_part,
              distribution = grid[i, "distribution"],
              inter = grid[i, "inter"]
            )
        },
        error = function(err) {
          print("error")
          model <- list(performance = "error")
        }
      )

      l <- list(cbind(grid[i, c("comb_id", "distribution", "inter")], model[, "performance"]))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)

    hyper_combinations <- lapply(hyper_combinations, function(x) dplyr::bind_rows(x)) %>%
      dplyr::bind_rows()

    if ("performance" %in% names(hyper_combinations)) {
      hyper_combinations <- hyper_combinations %>%
        dplyr::select(-performance)
    }

    hyper_combinations <- hyper_combinations %>%
      stats::na.omit()

    row.names(hyper_combinations) <- NULL

    ranked_combinations <- model_selection(hyper_combinations, metrics)

    # fit final model

    choosen_family <- ranked_combinations[[1]][1, "distribution"]
    full_data <- data
    if (families_bank[which(families_bank$distribution == choosen_family), "discrete"] == 1) {
      full_data[, "ind_ha"] <- round(full_data[, "ind_ha"])
    }

    message("\nFitting the best model...")
    final_model <-
      fit_abund_gam(
        data = full_data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        distribution = choosen_family,
        inter = ranked_combinations[[1]][1, "inter"]
      )

    message(
      "The best model was achieved with:",
      "\n family = ",
      choosen_family,
      "\n inter = ",
      ranked_combinations[[1]][1, "inter"]
    )

    final_list <- c(final_model, ranked_combinations)

    return(final_list)
  }

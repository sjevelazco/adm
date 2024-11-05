#' Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param grid tibble or data.frame. A dataframe with "family_call", "inter" as columns and its values combinations as rows.
#' @param metrics character. Vector with one or more metrics from c("corr_spear","corr_pear","mae","pdisp","inter","slope").
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
tune_abund_gam <-
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
    . <- discrete <- i <- performance <- NULL

    if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }

    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f)) %>%
      t() %>%
      as.vector()

    families_bank <-
      system.file("external/families_bank.txt", package = "adm") %>%
      utils::read.delim(., header = TRUE, quote = "\t") # families_bank

    # making grid
    if (is.null(grid)) {
      message("Families not provided. Picking recommended families to fit data.")
      family_call <- family_selector(data, response)$family_call
      inter <- "automatic"
      grid <- expand.grid(family_call = family_call, inter = inter)
      grid <- dplyr::left_join(grid, families_bank, by = "family_call") %>%
        dplyr::select(family_call, discrete, inter)
    } else if (is.data.frame(grid) & all(names(grid) %in% c("family_call", "inter")) & all(grid$family_call %in% families_bank$family_call)) {
      message("Testing with provided grid.")
      grid <- dplyr::left_join(grid, families_bank, by = "family_call") %>%
        dplyr::select(family_call, discrete, inter)
    } else {
      stop("Grid names expected to be 'family_call' and 'inter'.")
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

    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_gam", "adm_eval"), .packages = c("dplyr", "gamlss")) %dopar% {
      data_fam <- data
      if (grid[i, "discrete"] == 1) {
        data_fam[, response] <- round(data[, response])
      }

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
              distribution = grid[i, "family_call"],
              inter = grid[i, "inter"]
            )
        },
        error = function(err) {
          print("error")
          model <- list(performance = "error")
        }
      )

      l <- list(cbind(grid[i, c("comb_id", "family_call", "inter")], model[,"performance"]))
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

    choosen_family <- ranked_combinations[[1]][1, "family_call"]
    full_data <- data
    if (families_bank[which(families_bank$family_call == choosen_family), "discrete"] == 1) {
      full_data[, "ind_ha"] <- round(full_data[, "ind_ha"])
    }

    message("Fitting the best model...")
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
      "The best model was a GAM with:",
      "\n family = ",
      choosen_family,
      "\n inter = ",
      ranked_combinations[[1]][1, "inter"]
    )

    final_list <- c(final_model, ranked_combinations)

    return(final_list)
  }

#' Fit and validate Generalized Boosted Regression models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#' @param distribution character. A string specifying the distribution to be used. See gbm::gbm documentation for details.
#' @param n.trees integer. The total number of trees to fit.
#' @param interaction.depth integer. The maximum depth of each tree. Default 5
#' @param n.minobsinnode integer. The minimum number of observations in the terminal nodes of the trees. Default 5
#' @param shrinkage numeric. The learning rate of the algorithm. Default 0.1
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom gbm gbm predict.gbm
#' @importFrom stats formula sd
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gbm" class object from gbm package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each replica and partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#' 
#' # Database with species abundance and x and y coordinates
#' data("sppabund")
#' 
#' # Extract data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>% 
#'   dplyr::select(-.part2, -.part3)
#' 
#' # Explore reponse variables
#' some_sp$ind_ha %>% range()
#' some_sp$ind_ha %>% hist()
#' 
#' # Here we balance number of absences
#' some_sp <- 
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio=0.2)
#' 
#' # Fit a GBM model
#' mgbm <- fit_abund_gbm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12","elevation","sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "gaussian",
#'   n.trees = 100,
#'   interaction.depth = 5,
#'   n.minobsinnode = 5,
#'   shrinkage = 0.1,
#'   predict_part = TRUE
#' )
#' 
#' mgbm
#' }
fit_abund_gbm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           distribution,
           n.trees = 100,
           interaction.depth = 5,
           n.minobsinnode = 5,
           shrinkage = 0.1,
           verbose = TRUE) {
    . <- mae <- pdisp <- NULL

    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- dplyr::bind_rows(c(c = predictors))
    }


    # Adequate database
    data <- adapt_df(
      data = data,
      response = response,
      predictors = predictors,
      predictors_f = predictors_f,
      partition = partition
    )


    # ---- Formula ----
    if (is.null(fit_formula)) {
      formula1 <- stats::formula(paste(response, "~", paste(c(
        predictors,
        predictors_f
      ), collapse = " + ")))
    } else {
      formula1 <- fit_formula
    }

    if (verbose) {
      message(
        "Formula used for model fitting:\n",
        Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
        "\n"
      )
    }

    # Fit models
    np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))

    part_pred_list <- list()
    eval_partial_list <- list()

    for (h in 1:np) {
      if (verbose) {
        message("Replica number: ", h, "/", np)
      }
      # out <- pre_tr_te(data, p_names, h)

      folds <- data %>%
        dplyr::pull(p_names[h]) %>%
        unique() %>%
        sort()

      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()

      for (j in 1:length(folds)) {
        if (verbose) {
          message("-- Partition number ", j, "/", length(folds))
        }
        train_set <- data[data[, p_names[h]] != folds[j], ]
        test_set <- data[data[, p_names[h]] == folds[j], ]

        set.seed(13)
        model <- gbm::gbm(
          formula = formula1,
          data = train_set,
          distribution = distribution,
          n.trees = n.trees,
          interaction.depth = interaction.depth,
          n.minobsinnode = n.minobsinnode,
          shrinkage = shrinkage,
          bag.fraction = 0.9
        )

        pred <- gbm::predict.gbm(model, test_set, type = "response")
        observed <- dplyr::pull(test_set, response)
        eval_partial[[j]] <- dplyr::tibble(
          model = "gbm",
          adm_eval(obs = observed, pred = pred)
        )

        if (predict_part) {
          part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
        }
      }

      # Create final database with parameter performance
      names(eval_partial) <- 1:length(folds)
      eval_partial <-
        eval_partial[sapply(eval_partial, function(x) !is.null(dim(x)))] %>%
        dplyr::bind_rows(., .id = "partition")
      eval_partial_list[[h]] <- eval_partial

      if (predict_part) {
        names(part_pred) <- 1:length(folds)
        part_pred <-
          part_pred[sapply(part_pred, function(x) !is.null(dim(x)))] %>%
          dplyr::bind_rows(., .id = "partition")
        part_pred_list[[h]] <- part_pred
      }
    }


    # fit final model with all data
    set.seed(13)
    full_model <- gbm::gbm(
      formula = formula1,
      data = data,
      distribution = distribution,
      n.trees = n.trees,
      interaction.depth = interaction.depth,
      n.minobsinnode = n.minobsinnode,
      shrinkage = shrinkage,
      bag.fraction = 0.9
    )


    # bind predicted evaluation
    eval_partial <- eval_partial_list %>%
      dplyr::bind_rows(.id = "replica") %>%
      dplyr::as_tibble()

    # bind predicted partition
    if (predict_part) {
      part_pred <- part_pred_list %>%
        dplyr::bind_rows(.id = "replica")
    } else {
      part_pred <- NULL
    }

    # Summarize performance
    eval_final <- eval_partial %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(dplyr::across(mae:pdisp, list(
        mean = mean,
        sd = stats::sd
      )), .groups = "drop")

    variables <- dplyr::bind_cols(
      data.frame(
        model = "gbm",
        response = response
      ),
      variables
    ) %>% as_tibble()

    # Final object
    data_list <- list(
      model = full_model,
      predictors = variables,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )
    return(data_list)
  }

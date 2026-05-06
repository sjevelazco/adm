#' Fit and validate Generalized Additive Models
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
#' @param inter integer. Number of knots in x-axis. Default "automatic"
#' @param distribution character. A string specifying the distribution to be used. See \link[gamlss.dist]{gamlss.family} documentation for details. Use distribution = gamlss.dist::NO(). Default NULL
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#' @param control_gamlss function. control parameters of the outer iterations algorithm in gamlss
#' See \link[gamlss]{gamlss.control} documentation for details. Default gamlss.control()
#'
#' @importFrom dplyr bind_rows bind_cols pull tibble as_tibble group_by summarise across
#' @importFrom gamlss gamlss pb predictAll gamlss.control
#' @importFrom stats formula sd
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gamlss" class object from gamlss package. This object can be used for predicting.
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
#' require(gamlss)
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
#'   balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)
#'
#' # Explore different family distributions
#' family_selector(data = some_sp, response = "ind_ha") %>% tail()
#'
#' # Fit a GAM model
#' mgam <- fit_abund_gam(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   sigma_formula = ~ elevation + bio3 + bio12,
#'   predictors_f = NULL,
#'   partition = ".part",
#'   distribution = gamlss.dist::ZAGA()
#' )
#'
#' mgam
#' }
fit_abund_gam <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           sigma_formula = ~1,
           nu_formula = ~1,
           tau_formula = ~1,
           partition,
           hold_out_set = NULL,
           predict_part = FALSE,
           distribution = NULL,
           inter = "automatic",
           verbose = TRUE,
           control_gamlss = gamlss::gamlss.control(trace = FALSE)) {
    . <- mae <- pdisp <- NULL

    if (is.null(distribution)) {
      stop("'distribution' argument was not used, a distribution must be specifyied")
    }

    # Adequate database
    data <- adapt_df(
      data = data,
      response = response,
      predictors = predictors,
      predictors_f = predictors_f,
      partition = partition
    )

    # Adequate hold-out set
    hold_out_set <- check_adapt_holdout_set(
      hold_out_set,
      predictors,
      predictors_f,
      response
    )
    hold_out_evaluation <- !is.null(hold_out_set)

    # Variables
    variables <- get_variables(predictors, predictors_f)

    # Formula
    if (is.null(fit_formula)) {
      if (inter == "automatic") {
        formula1 <-
          paste(c(
            paste("pb(", predictors, ")", collapse = " + ", sep = ""),
            predictors_f
          ), collapse = " + ")
        formula1 <- stats::formula(paste(
          response, "~", formula1
        ))
      } else {
        formula1 <-
          paste(c(
            paste("pb(", predictors, paste0(", inter = ", inter, ")"), collapse = " + ", sep = ""),
            predictors_f
          ), collapse = " + ")
        formula1 <- stats::formula(paste(
          response, "~", formula1
        ))
      }
    } else {
      formula1 <- fit_formula
    }

    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
      "\n"
    )


    # Fit models
    if (is.null(partition) || !any(nzchar(partition, keepNA = FALSE))) {
      set.seed(13)
      suppressWarnings(
        full_model <- gamlss::gamlss(
          formula = formula1,
          family = distribution,
          data = data,
          sigma.formula = sigma_formula,
          nu.formula = nu_formula,
          tau.formula = tau_formula,
          control = control_gamlss,
          trace = FALSE
        )
      )
      result <- list(
        model = full_model
      )
      return(result)
    } else {
      np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
      p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))

      # part_pred_list <- list()
      # eval_partial_list <- list()

      replica_training_lists <- init_training_lists("replica")

      family <- distribution
      for (h in 1:np) {
        if (verbose) {
          message("Replica number: ", h, "/", np)
        }

        folds <- data %>%
          dplyr::pull(p_names[h]) %>%
          unique() %>%
          sort()

        fold_training_lists <- init_training_lists("fold")

        # eval_partial <- list()
        # pred_test <- list()
        # part_pred <- list()

        for (j in 1:length(folds)) {
          if (verbose) {
            message("-- Partition number ", j, "/", length(folds))
          }

          train_set <- data[data[, p_names[h]] != folds[j], ]
          test_set <- data[data[, p_names[h]] == folds[j], ]

          set.seed(13)
          model <- gamlss::gamlss(
            formula = formula1,
            family = family,
            data = train_set,
            sigma.formula = sigma_formula,
            nu.formula = nu_formula,
            tau.formula = tau_formula,
            control = control_gamlss,
            trace = FALSE
          )

          pred <- gamlss::predictAll(model, newdata = test_set, data = train_set, type = "response")[[1]]
          observed <- dplyr::pull(test_set, response)

          if (hold_out_evaluation) {
            pred_ho <-
              suppressMessages(stats::predict(model, newdata = hold_out_set[, c(predictors, predictors_f)], type = "response"))
            observed_ho <- hold_out_set[, response]
          } else {
            pred_ho <- observed_ho <- NULL
          }

          fold_training_lists <- fold_perf_register(
            "gam", folds, j,
            fold_training_lists,
            predict_part,
            hold_out_evaluation,
            pred, pred_ho,
            observed, observed_ho
          )
        }

        # Create final database with parameter performance
        replica_training_lists <- replica_perf_register(
          replica_training_lists, fold_training_lists,
          folds, h, predict_part, hold_out_evaluation
        )
      }


      # fit final model with all data
      set.seed(13)
      full_model <- gamlss::gamlss(
        formula = formula1,
        family = family,
        data = data,
        sigma.formula = sigma_formula,
        nu.formula = nu_formula,
        tau.formula = tau_formula,
        control = control_gamlss,
        trace = FALSE
      )

      # evaluate full model with hold-out set
      if (hold_out_evaluation) {
        pred <-
          suppressMessages(predict(full_model, newdata = hold_out_set[, c(predictors, predictors_f)], type = "response"))
        observed <- hold_out_set[, response]

        hold_out_perf <- adm_eval(obs = observed, pred = pred)
      } else {
        hold_out_perf <- NULL
      }

      # Construct the standard final list to be returned
      data_list <- wrap_final_list(
        "gam",
        full_model,
        variables,
        response,
        replica_training_lists,
        hold_out_evaluation,
        hold_out_perf,
        predict_part,
        get_metadata(
          "gam",
          list(
            formula = formula1,
            sigma.formula = sigma_formula,
            nu.formula = nu_formula,
            tau.formula = tau_formula
          )
        )
      )

      return(data_list)
    }
  }

#' Fit and validate Support Vector Machine models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative
#' predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor
#' variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param sigma numeric. Inverse kernel width for the Radial Basis kernel function "rbfdot". Default "automatic".
#' @param C numeric. Cost of constraints violation, the 'C'-constant of the regularization term in the Lagrange formulation. Default 1
#' @param predict_part logical. Save predicted abundance for testing data. Default = FALSE
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "ksvm" class object from kernlab package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#' @export
#'
#' @examples
fit_abund_svm <-
  function(data,
           response,
           predictors,
           predictors_f,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           sigma = "automatic",
           C = 1) {
    # Formula
    if (is.null(fit_formula)) {
      formula1 <- stats::formula(paste(response, "~", paste(c(
        predictors,
        predictors_f
      ), collapse = " + ")))
    } else {
      formula1 <- fit_formula
    }

    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()

    eval_partial <- list()
    part_pred <- list()

    if (sigma == "automatic") {
      kpar_ <- "automatic"
    } else {
      kpar_ <- list(sigma = sigma)
    }

    for (j in 1:length(folds)) {
      message("-- Evaluating with fold ", j, "/", length(folds))

      train_set <- data[data[, partition] != folds[j], ]
      test_set <- data[data[, partition] == folds[j], ]

      test_rf <- kernlab::ksvm(
        formula1,
        data = train_set,
        type = "eps-svr", # TODO check nu-svr
        kernel = "rbfdot",
        kpar = kpar_,
        C = C
      )

      pred <- kernlab::predict(test_rf, test_set, type = "response")
      observed <- dplyr::pull(test_set, response)
      eval_partial[[j]] <- dplyr::tibble(
        model = "svm",
        adm_eval(obs = observed, pred = pred)
      )

      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }

    # fit final model with all data
    model <- kernlab::ksvm(
      formula1,
      data = data,
      type = "eps-svr",
      kernel = "rbfdot",
      kpar = kpar_,
      C = C
    )


    # bind predicted evaluation
    eval_partial <- eval_partial %>%
      dplyr::bind_rows() %>%
      dplyr::as_tibble()

    # bind predicted partition
    if (predict_part) {
      part_pred <- part_pred %>%
        dplyr::bind_rows() %>%
        dplyr::as_tibble()
    } else {
      part_pred <- NULL
    }

    # Sumarize performance
    eval_final <- eval_partial %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(dplyr::across(corr_spear:pdisp, list(
        mean = mean,
        sd = stats::sd
      )), .groups = "drop")

    # Final object
    data_list <- list(
      model = model,
      performance = eval_final,
      performance_part = eval_partial,
      predicted_part = part_pred
    )
    return(data_list)
  }

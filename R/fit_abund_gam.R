#' fit_abund_gam
#'
#' @param data
#' @param response
#' @param predictors
#' @param predictors_f
#' @param fit_formula
#' @param partition
#' @param predict_part
#' @param k
#' @param family
#' @param method
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom mgcv gam
#' @importFrom stats formula predict sd
#'
#' @return
#' @export
#'
#' @examples
fit_abund_gam <-
  function(data,
           response,
           predictors,
           predictors_f,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           family,
           inter) {
    # Variables
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    # Formula
    if (is.null(fit_formula)) {
      formula1 <-
        paste(c(
          paste("pb(", predictors, paste0(", inter = ", inter, ")"), collapse = " + ", sep = ""),
          predictors_f
        ), collapse = " + ")
      formula1 <- stats::formula(paste(
        response, "~", formula1
      ))
    } else {
      formula1 <- fit_formula
    }

    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()

    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      message("-- Evaluating with fold ", j, "/", length(folds))

      train_set <- data[data[, partition] != folds[j], ]
      test_set <- data[data[, partition] == folds[j], ]

      model <- gamlss::gamlss(
        formula = formula1,
        family = family,
        data = train_set
      )

      pred <- predict(model, newdata = test_set, data = train_set, type = "response")
      observed <- dplyr::pull(test_set, response)
      eval_partial[[j]] <- dplyr::tibble(
        model = "gam",
        adm_eval(obs = observed, pred = pred)
      )

      if (predict_part) {
        part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
      }
    }
    
    # fit final model with all data
    full_model <- gamlss::gamlss(
      formula = formula1,
      family = family,
      data = data
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

    # Summarize performance
    eval_final <- eval_partial %>%
      dplyr::group_by(model) %>%
      dplyr::summarise(dplyr::across(corr_spear:pdisp, list(
        mean = mean,
        sd = stats::sd
      )), .groups = "drop")

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

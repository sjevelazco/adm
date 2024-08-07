#' Fit and validate Generalized Linear Models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param family character. A string specifying the distribution to be used. See gamlss::gamlss documentation for details.
#' @param poly interger >= 2. If used with values >= 2 model will use polynomials for those continuous variables (i.e. used in predictors argument). Default is 0.
#' @param inter_order interger >= 0. The interaction order between explanatory variables. Default is 0.
#' @param verbose logical. If FALSE, disables all console messages. Default TRUE
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom gamlss gamlss
#' @importFrom stats formula sd
#'
#' @return
#'
#' A list object with:
#' \itemize{
#' \item model: A "gamlss" class object from gamlss package. This object can be used for predicting.
#' \item predictors: A tibble with quantitative (c column names) and qualitative (f column names) variables use for modeling.
#' \item performance: Averaged performance metrics (see \code{\link{adm_eval}}).
#' \item performance_part: Performance metrics for each partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
#'
#' @export
#'
#' @examples
fit_abund_glm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           family,
           poly = 0,
           inter_order = 0,
           verbose = TRUE) {
    # Variables
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

    # TODO Formula
    if (is.null(fit_formula)) {
      if (poly >= 2) {
        forpoly <- lapply(2:poly, function(x) {
          paste("I(", predictors, "^", x, ")",
            sep = "", collapse = " + "
          )
        }) %>% paste(collapse = " + ")
        formula1 <- paste(c(predictors, predictors_f), collapse = " + ") %>%
          paste(., forpoly, sep = " + ")
      } else {
        formula1 <-
          paste(c(predictors, predictors_f), collapse = " + ")
      }

      if (inter_order > 0) {
        forinter <- c(predictors, predictors_f)
        if (inter_order > length(forinter)) {
          stop("value of inter_order is higher than number of predicors ", "(", length(forinter), ")")
        }
        forinter_l <- list()

        for (i in 1:inter_order) {
          forinter_l[[i]] <- do.call(
            "expand.grid",
            c(lapply(1:(i + 1), function(x) {
              forinter
            }), stringsAsFactors = FALSE)
          )
          forinter_l[[i]] <- apply(forinter_l[[i]], 1, function(x) {
            x <- unique(sort(x))
            if (length(x) > i) {
              paste(x, collapse = ":")
            }
          }) %>%
            unlist() %>%
            unique()
        }
        forinter <- sapply(forinter_l, paste, collapse = " + ")
        forinter <- do.call("paste", c(as.list(forinter), sep = " + "))
      }

      if (exists("forinter")) {
        formula1 <- paste(formula1, forinter, sep = " + ")
        formula1 <- stats::formula(paste(
          response, "~", formula1
        ))
      } else {
        formula1 <- stats::formula(paste(
          response, "~", formula1
        ))
      }
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

    folds <- data %>%
      dplyr::pull(partition) %>%
      unique() %>%
      sort()

    eval_partial <- list()
    part_pred <- list()
    for (j in 1:length(folds)) {
      if (verbose) {
        message("-- Evaluating with fold ", j, "/", length(folds))
      }


      # if (family=="poisson"){
      #   data[, response] <- round(data[, response], 0)
      #   data[, response] <- as.integer(data[[response]])
      # }

      train_set <- data[data[, partition] != folds[j], ]
      test_set <- data[data[, partition] == folds[j], ]

      model <- gamlss::gamlss(
        formula = formula1,
        family = family,
        data = train_set,
        trace = FALSE
      )

      pred <- predict(model, newdata = test_set, data = train_set, type = "response")
      observed <- dplyr::pull(test_set, response)
      eval_partial[[j]] <- dplyr::tibble(
        model = "glm",
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
      data = data,
      trace = FALSE
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

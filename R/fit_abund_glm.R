#' Fit and validate Generalized Linear Models
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
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param distribution character. A string specifying the distribution to be used. See \link[gamlss.dist]{gamlss.family} documentation for details. Use distribution = gamlss.dist::NO(). Default NULL
#' @param poly integer >= 2. If used with values >= 2 model will use polynomials for those continuous variables (i.e. used in predictors argument). Default is 0.
#' @param inter_order integer >= 0. The interaction order between explanatory variables. Default is 0.
#' @param control_gamlss function. control parameters of the outer iterations algorithm in gamlss
#' See \link[gamlss]{gamlss.control} documentation for details. Default gamlss.control()
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
#' \item performance_part: Performance metrics for each replica and partition.
#' \item predicted_part: Observed and predicted abundance for each test partition.
#' }
#'
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
#' # Fit a GLM model
#' glm_1 <- fit_abund_glm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "ZAGA",
#'   poly = 0,
#'   inter_order = 0,
#'   predict_part = TRUE
#' )
#'
#' glm_1
#'
#' # Using second order polynomials and first order interaction terms
#' glm_2 <- fit_abund_glm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "ZAGA",
#'   poly = 2,
#'   inter_order = 1,
#'   predict_part = TRUE
#' )
#'
#' glm_2
#'
#' # Using third order polynomials and second order interaction terms
#' glm_3 <- fit_abund_glm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "ZAGA",
#'   poly = 3,
#'   inter_order = 2,
#'   predict_part = TRUE
#' )
#'
#' glm_3
#'
#' # Setting formulas for different distribution parameters
#' glm_4 <- fit_abund_glm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "ZAGA",
#'   fit_formula = ind_ha ~ bio12 + elevation + sand + eco,
#'   sigma_formula = ind_ha ~ bio12 + elevation + sand,
#'   poly = 0,
#'   inter_order = 0,
#'   predict_part = TRUE
#' )
#'
#' glm_4
#' }
fit_abund_glm <-
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
           distribution = NULL,
           poly = 0,
           inter_order = 0,
           control_gamlss = gamlss::gamlss.control(trace = FALSE),
           verbose = TRUE) {
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


    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- dplyr::bind_rows(c(c = predictors))
    }


    # Formula
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

    # Fit models
    np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))

    part_pred_list <- list()
    eval_partial_list <- list()
    family <- distribution

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

        pred <- predict(model, newdata = test_set, data = train_set, type = "response")
        observed <- dplyr::pull(test_set, response)
        eval_partial[[j]] <- dplyr::tibble(
          model = "glm",
          adm_eval(obs = observed, pred = pred)
        )

        if (predict_part) {
          part_pred[[j]] <- dplyr::tibble(partition = folds[j], observed, predicted = pred)
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
        model = "glm",
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

    # Standardize output list
    for (i in 2:length(data_list)) {
      if (!class(data_list[[i]])[1] == "tbl_df") {
        data_list[[i]] <- dplyr::as_tibble(data_list[[i]])
      }
    }

    return(data_list)
  }

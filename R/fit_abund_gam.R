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
#' See \link[gamlss.contro]{gamlss} documentation for details. Default gamlss.control()
#'
#' @importFrom dplyr bind_rows pull tibble as_tibble group_by summarise across
#' @importFrom gamlss gamlss pb predictAll
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
           predict_part = FALSE,
           distribution = NULL,
           inter = "automatic",
           verbose = TRUE,
           control_gamlss = gamlss.control(trace = FALSE)
           ) {
    . <- mae <- pdisp <-  NULL
    
    if (is.null(distribution)) {
      stop("'distribution' argument was not used, a distribution must be specifyied")
    }
    
    # Adequate database
    data <- adapt_df(data = data,
                     response = response,
                     predictors = predictors,
                     predictors_f = predictors_f, 
                     partition = partition)
    
    # Variables
    if (!is.null(predictors_f)) {
      variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
    } else {
      variables <- dplyr::bind_rows(c(c = predictors))
    }
    
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
    np <- ncol(data %>% dplyr::select(dplyr::starts_with(partition)))
    p_names <- names(data %>% dplyr::select(dplyr::starts_with(partition)))
    
    part_pred_list <- list()
    eval_partial_list <- list()
    
    for (h in 1:np) {
      if (verbose) {
        message("Replica number: ", h, "/", np)
      }
      # out <- pre_tr_te(data, p_names, h)
      
      folds <- data %>% dplyr::pull(p_names[h]) %>% unique() %>% sort()
      
      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()
      family <- distribution
      
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
        eval_partial[[j]] <- dplyr::tibble(
          model = "gam",
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

    variables <- bind_cols(
      data.frame(
        model = "gam",
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

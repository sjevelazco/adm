#' Fit and validate Generalized Linear Models
#'
#' @param data tibble or data.frame. Database with response, predictors, and partition values
#' @param response character. Column name with species abundance.
#' @param predictors character. Vector with the column names of quantitative predictor variables (i.e. continuous variables). Usage predictors = c("temp", "precipt", "sand")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor variables (e.g. formula(abund ~ temp + precipt + sand + landform)). Note that the variables used here must be consistent with those used in response, predictors, and predictors_f arguments. Default NULL
#' @param partition character. Column name with training and validation partition groups.
#' @param predict_part logical. Save predicted abundance for testing data. Default is FALSE.
#' @param distribution character. A string specifying the distribution to be used. See \link[gamlss.dist]{gamlss.family} documentation for details. Use distribution = gamlss.dist::NO(). Default NULL
#' @param poly integer >= 2. If used with values >= 2 model will use polynomials for those continuous variables (i.e. used in predictors argument). Default is 0.
#' @param inter_order integer >= 0. The interaction order between explanatory variables. Default is 0.
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
fit_abund_glm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           distribution = NULL,
           poly = 0,
           inter_order = 0,
           verbose = TRUE) {
    
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
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))

   
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
      
      folds <- data %>% dplyr::pull(p_names[h]) %>% unique() %>% sort()
      
      eval_partial <- list()
      pred_test <- list()
      part_pred <- list()
      
      for (j in 1:length(folds)) {
        if (verbose) {
          message("-- Partition number ", j, "/", length(folds))
        }
        
        train_set <- data[data[, p_names[h]] != folds[j], ]
        test_set <- data[data[, p_names[h]] == folds[j], ]
        
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
    full_model <- gamlss::gamlss(
      formula = formula1,
      family = family,
      data = data,
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

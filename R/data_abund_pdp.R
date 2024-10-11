#' Calculate data to construct partial dependence plots
#'
#' @param model object returned by any fit_abund or tune_abund family functions
#' @param predictors character. Vector with two predictor name(s) to plot. If NULL all predictors will be plotted. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict continuous predictors. Default 50
#' @param resid logical. Calculate residuals based on training data. Default FALSE
#' @param training_data data.frame. Database with response and predictor values used
#' to fit a model. Default NULL. Required for GLM, GAM, DNN, NET, RAF, SVM models
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. When this argument is used, function will calculate partial dependence curves
#' distinguishing conditions used in training and projection conditions
#' (i.e., projection data present in projection area but not training). Default NULL
#'
#' @importFrom dplyr select tibble
#' @importFrom gbm predict.gbm
#' @importFrom kernlab predict
#' @importFrom stats na.omit
#' @importFrom terra minmax
#'
#' @return A list with two tibbles "pdpdata" and "resid".
#' \itemize{
#' \item pdpdata: has data to construct partial dependence plots, the first column includes values of the selected environmental
#' variable, the second column with predicted suitability, and the third
#'  column with range type, with two values Training and Projecting, referring to suitability
#'  calculated within and outside the range of training conditions. Third column is only returned
#'  if "projection_data" argument is used
#' \item resid: has data to plot residuals. The first column includes values of the selected environmental
#'  variable and the second column with predicted suitability.
#' }
#'
#' @export
#'
#' @examples
data_abund_pdp <-
  function(model,
           predictors,
           resolution = 50,
           resid = FALSE,
           training_data = NULL,
           invert_transform = NULL,
           response_name = NULL,
           projection_data = NULL) {
    self <- NULL

    # Extract training data
    # TODO
    # if (class(model)[1] == "gam") {
    #   x <- model$model[attr(model$terms, "term.labels")]
    # }

    # TODO
    # if (class(model)[1] == "graf") {
    #   x <- model$obsx
    #   x <- x[names(model$peak)]
    # }

    # TODO
    # if (class(model)[1] == "glm") {
    #   flt <- grepl("[I(]", attr(model$terms, "term.labels")) |
    #     grepl(":", attr(model$terms, "term.labels"))
    #   flt <- attr(model$terms, "term.labels")[!flt]
    #   x <- model$model[flt]
    # }

    if (class(model)[1] == "list") {
      if (all(c("model", "predictors", "performance", "performance_part", "predicted_part") %in% names(model))
      ) {
        variables <- model$predictors
        model <- model[[1]]
      }
    }

    if (any(class(model)[1] == c("gamlss", "luz_module_fitted"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for Generalized Linear Models (GLM), Generalized Additive Models (GAM) and Deep Neural Network (DNN) it is necessary to provide calibration data in 'training_data' argument"
        )
      }
      x <- training_data[, as.vector(variables[1, ])[2:ncol(variables)] %>% unlist()]
    }

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula", "ksvm", "gbm"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for Neural Networks (NET), Random Forest (RAF), Support Vector Machine (SVM) it is necessary to provide calibration data in 'training_data' argument"
        )
      }

      if (class(model)[1] == "ksvm") {
        x <- training_data[names(attr(model@terms, "dataClasses")[-1])]
      } else if (class(model)[1] == "gbm") {
        x <- training_data[, c(model$response.name, model$var.names)]
      } else {
        x <- training_data[names(attr(model$terms, "dataClasses")[-1])]
      }
    }

    x <- stats::na.omit(x)

    # Control average factor level
    fact <- sapply(x, is.factor)
    suit_c <- which(!fact) # TODO suit_c tem algo a ver com suitability?
    fact <- which(fact)
    suit_c <- data.frame(t(apply(x[suit_c], 2, mean)))
    # suit_c <- data.frame((x[1, ])) For residuals

    if (sum(fact) > 0) {
      for (i in 1:length(fact)) {
        ff <- sort(data.frame(unique(x[names(fact[i])]))[, 1])
        ff <- ff[as.integer(length(ff) / 2)]
        suit_c[names(fact)[i]] <- ff
      }
    }


    if (predictors %in% names(fact)) {
      rng <- sort(data.frame(unique(x[names(fact)]))[, predictors])
      suit_c <- suit_c %>% dplyr::select(-{{ predictors }})
    } else {
      if (is.null(projection_data)) {
        rng <- range(x[, predictors])
        rng <- seq(rng[1], rng[2], length.out = resolution)
      } else {
        # Range extrapolation
        rng <- terra::minmax(projection_data[[predictors]])
        rng <- seq(rng[1], rng[2], length.out = resolution)
      }
    }

    suit_c <- data.frame(rng, suit_c)
    suit_c[predictors] <- NULL
    names(suit_c)[1] <- predictors

    # Predict model

    if (class(model)[1] == "luz_module_fitted") {
      create_dataset <- torch::dataset(
        "dataset",
        initialize = function(df, response_variable = 0) {
          self$df <- df
        },
        .getitem = function(index) {
          x <- torch::torch_tensor(as.numeric(self$df[index, ]))
          list(x = x)
        },
        .length = function() {
          nrow(self$df)
        }
      )

      pred_dataset <- create_dataset(suit_c %>% dplyr::select(-variables[["response"]]))
      pred_dataset_x <- create_dataset(x %>% dplyr::select(-variables[["response"]]))

      suit_c <-
        data.frame(suit_c[1],
          Abundance = suppressMessages(
            stats::predict(
              model,
              newdata = pred_dataset,
              type = "response"
            ) %>%
              as.numeric()
          )
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
            Abundance = suppressMessages(
              stats::predict(
                model,
                newdata = pred_dataset_x,
                type = "response"
              ) %>%
                as.numeric()
            )
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "gamlss") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance = suppressMessages(
            predict(
              model,
              newdata = suit_c,
              data = training_data,
              type = "response"
            )
          )
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Abundance = suppressMessages(predict(model, newdata = x, data = training_data, type = "response")))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "gbm") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance = suppressMessages(gbm::predict.gbm(model, newdata = suit_c, type = "response"))
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors], Abundance = suppressMessages(gbm::predict.gbm(model, newdata = x, type = "response")))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }


    if (class(model)[1] == "nnet.formula") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance =
            suppressMessages(stats::predict(
              model,
              newdata = suit_c, type = "raw"
            ))
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
            Abundance =
              suppressMessages(stats::predict(
                model,
                type = "raw"
              ))
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "randomForest.formula") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance =
            suppressMessages(stats::predict(model, suit_c, type = "response"))
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
            Abundance =
              suppressMessages(stats::predict(model, type = "response"))
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance = kernlab::predict(model, suit_c, type = "response")[, 2]
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
            Abundance = kernlab::predict(model, x, type = "response")[, 2]
          )
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    # Category of training and projection data
    if (!is.null(projection_data)) {
      if (!predictors %in% names(fact)) {
        result$pdpdata$Type <-
          ifelse(suit_c[, 1] >= min(x[, predictors]) &
            suit_c[, 1] <= max(x[, predictors]),
          "Training",
          "Projection"
          )
      }
    }

    result <- lapply(result, function(x) if (is.data.frame(x)) dplyr::tibble(x))

    if (!is.null(invert_transform)) {
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) {
          next
        } else {
          result[[i]] <- result[[i]] %>% dplyr::mutate(
            Abundance = adm_transform(result[[i]],
              "Abundance",
              invert_transform[["method"]],
              inverse = TRUE,
              t_terms = c(
                invert_transform[["a"]] %>% as.numeric(),
                invert_transform[["b"]] %>% as.numeric()
              )
            ) %>% dplyr::pull(Abundance_inverted)
          )
        }
      }
    }

    if (!is.null(response_name)) {
      for (i in 1:length(result)) {
        if (is.null(result[[i]])) {
          next
        } else {
          idx <- which(names(result[[i]]) == "Abundance")
          names(result[[i]])[[idx]] <- response_name
        }
      }
    }

    return(result)
  }

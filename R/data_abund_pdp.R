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
#' @param invert_transform vector. A vector containing method and terms to invert transformation of response variable. Useful for those cases that the response variable was transformed with one of the method in \code{\link{adm_transform}}. Usage:
#' \itemize{
#' \item For "01": invert_transform = c(method = "01", a = min(x), b = max(x))
#' \item For "zscore": invert_transform = c(method = "zscore", a = mean(x), b = sd(x))
#' \item For "log" and "log1: not needed.
#' \item Can't invert "round" transformations.
#' \item Default NULL
#' }
#' @param response_name character. Name of the response variable. Default NULL
#' @param sample_size vector. For CNN only. A vector containing the dimensions, in pixels, of raster samples. See cnn_make_samples beforehand. Default c(11,11)
#' @param training_raster a terra SpatRaster object. For CNN only. A raster containing the predictor variables used in tune_abund_cnn or fit_abund_cnn.
#' @param x_coord character. For CNN only. The name of the column containing longitude information for each observation.
#' @param y_coord character. For CNN only. The name of the column containing latitude information for each observation.
#'
#' @importFrom dplyr select all_of bind_cols tibble mutate pull
#' @importFrom gbm predict.gbm
#' @importFrom kernlab predict
#' @importFrom stats na.omit model.matrix
#' @importFrom terra minmax spatSample
#' @importFrom torch torch_manual_seed dataset torch_tensor
#' @importFrom torchvision transform_to_tensor
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
#' @seealso \code{\link{data_abund_bpdp}}, \code{\link{p_abund_pdp}}, \code{\link{p_abund_bpdp}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' # Load  data
#' envar <- system.file("external/envar.tif", package = "adm") %>%
#'   rast()
#'
#' data("sppabund")
#' some_sp <- sppabund %>%
#'   filter(species == "Species one")
#'
#' # Fit some models
#' mglm <- fit_abund_glm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("bio12", "elevation", "sand"),
#'   predictors_f = c("eco"),
#'   partition = ".part",
#'   distribution = "ZAIG",
#'   poly = 3,
#'   inter_order = 0,
#'   predict_part = TRUE
#' )
#'
#' # Prepare data for Partial Dependence Plots
#' pdp_data <- data_abund_pdp(
#'   model = mglm,
#'   predictors = "bio12",
#'   resolution = 25,
#'   resid = TRUE,
#'   training_data = some_sp,
#'   response_name = "Abundance",
#'   projection_data = envar
#' )
#'
#' pdp_data
#' }
data_abund_pdp <-
  function(model,
           predictors,
           resolution = 50,
           resid = FALSE,
           training_data = NULL,
           invert_transform = NULL,
           response_name = NULL,
           projection_data = NULL,
           sample_size = c(11, 11),
           training_raster = NULL,
           x_coord = NULL,
           y_coord = NULL) {
    self <- Abundance_inverted <- y <- NULL

    if (length(predictors) > 1) {
      stop("Please provide only one predictor.")
    }

    if (class(model)[1] == "list") {
      if (all(c("model", "predictors", "performance", "performance_part", "predicted_part") %in% names(model))
      ) {
        variables <- model$predictors
        model <- model[[1]]
      }
    } else {
      stop('Please, use tune_abund_ or fit_abund_ output list in "model" argument.')
    }

    # Check if the required parameters for cnn
    if (class(model)[1] == "luz_module_fitted" & variables[["model"]] == "cnn") {
      if (is.null(sample_size)) {
        stop("sample_size is needed. Use the same as in tune_abund_cnn or fit_abund_cnn")
      } else if (is.null(training_raster)) {
        stop("training_raster is needed. Use the same as in tune_abund_cnn or fit_abund_cnn")
      } else if (is.null(x_coord) | is.null(y_coord)) {
        stop("x_coord and y_coord are needed. Use the x and y arguments of tune_abund_cnn or fit_abund_cnn")
      }
    }

    if (any(class(model)[1] == c("gamlss", "luz_module_fitted", "xgb.Booster"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for Generalized Linear Models (GLM), Generalized Additive Models (GAM) and Deep Neural Network (DNN) it is necessary to provide calibration data in 'training_data' argument"
        )
      }
      if (variables[["model"]] == "cnn") {
        x <- training_data[, c(as.vector(variables[1, ])[2:ncol(variables)] %>% unlist(), x_coord, y_coord)]
      } else {
        x <- training_data[, as.vector(variables[1, ])[2:ncol(variables)] %>% unlist()]
      }
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
    suit_c <- which(!fact)
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

    # Make cnn samples
    if (variables[["model"]] == "cnn") {
      suit_c <- suit_c %>%
        dplyr::select(-dplyr::all_of(x_coord), -dplyr::all_of(y_coord))

      training_raster <- training_raster[[names(suit_c %>% dplyr::select(-variables[["response"]]))]]

      for (n in names(training_raster)) {
        if (n != predictors) {
          training_raster[[n]] <- suit_c[[n]] %>% unique()
        }
      }

      randompoints <- terra::spatSample(training_raster, resolution, xy = T) %>%
        dplyr::select(x, y) %>%
        dplyr::bind_cols(suit_c)

      random_samples <- cnn_make_samples(
        data = randompoints,
        x = "x",
        y = "y",
        response = variables[["response"]],
        size = cnn_get_crop_size(sample_size),
        raster = training_raster,
        raster_padding = TRUE,
        padding_method = "zero"
      )

      for (i in 1:length(random_samples$predictors)) {
        random_samples$predictors[[i]][, , predictors] <- suit_c[i, predictors]
      }

      x_samples <- cnn_make_samples(
        data = x,
        x = x_coord,
        y = y_coord,
        response = variables[["response"]],
        size = cnn_get_crop_size(sample_size),
        raster = training_raster
      )
    }

    # Predict model

    #### cnn ####
    if (class(model)[1] == "luz_module_fitted" & variables[["model"]] == "cnn") {
      model$model$eval()
      torch::torch_manual_seed(13)

      create_dataset <- torch::dataset(
        "dataset",
        initialize = function(data_list) {
          self$predictors <- data_list$predictors
        },
        .getitem = function(index) {
          x <- torchvision::transform_to_tensor(self$predictors[[index]])
          list(x = x)
        },
        .length = function() {
          length(self$predictors)
        }
      )

      pred_dataset <- create_dataset(random_samples)
      pred_dataset_x <- create_dataset(x_samples)

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

    #### dnn ####
    if (class(model)[1] == "luz_module_fitted" & variables[["model"]] == "dnn") {
      model$model$eval()
      torch::torch_manual_seed(13)

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

    #### gam and glm ####
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

    #### xgb ####
    if (class(model)[1] == "xgb.Booster") {
      pred_matrix <- list(
        data = stats::model.matrix(~ . - 1, data = suit_c[, model$feature_names])
      )

      suit_c <-
        data.frame(suit_c[1],
          Abundance = suppressMessages(stats::predict(model, newdata = pred_matrix$data, type = "response"))
        )

      if (resid) {
        pred_matrix_r <- list(
          data = stats::model.matrix(~ . - 1, data = x[, model$feature_names])
        )

        suit_r <-
          data.frame(x[predictors], Abundance = suppressMessages(stats::predict(model, newdata = pred_matrix_r$data, type = "response")))
        result <- list("pdpdata" = suit_c, "resid" = suit_r)
      } else {
        result <- list("pdpdata" = suit_c, "resid" = NA)
      }
    }

    #### gbm ####
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

    #### net ####
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

    #### raf ####
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

    #### svm ####
    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1],
          Abundance = kernlab::predict(model, suit_c, type = "response")
        )
      if (resid) {
        suit_r <-
          data.frame(x[predictors],
            Abundance = kernlab::predict(model, x, type = "response")
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

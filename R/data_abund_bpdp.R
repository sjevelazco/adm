#' Calculate data to construct bivariate partial dependence plots
#'
#' @description Calculate data to construct bivariate partial dependence for two predictor set
#'
#' @param model object returned by any fit_abund or tune_abund family functions
#' @param predictors character. Vector with two predictor name(s) to plot. If NULL all predictors will be plotted. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict continuous predictors. Default 50
#' @param training_data data.frame or tibble. Database with response (0,1) and predictor values used to fit a model. Default NULL
#' @param training_boundaries character. Plot training conditions boundaries based on training data (i.e., presences, presences and absences, etc). If training_boundaries = "convexh", function will delimit training environmental region based on a convex-hull. If training_boundaries = "rectangle", function will delimit training environmental region based on four straight lines. If used any methods it is necessary provide data in training_data argument. If NULL all predictors will be used. Default NULL.
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model projection. Default NULL
#' @param invert_transform logical. Invert transformation of response variable. Useful for those cases that the response variable was transformed with one of the method in \code{\link{adm_transform}}. Default NULL
#' @param response_name character. Name of the response variable. Default NULL
#'
#' @importFrom dplyr as_tibble select
#' @importFrom gbm predict.gbm
#' @importFrom grDevices chull
#' @importFrom kernlab predict
#' @importFrom stats na.omit
#' @importFrom terra minmax
#' @importFrom torch dataset torch_tensor
#'
#' @seealso \code{\link{data_abund_pdp}}, \code{\link{p_abund_pdp}}, \code{\link{p_abund_bpdp}}
#'
#' @return A list with two tibbles "pdpdata" and "resid".
#' \itemize{
#' \item pdpdata: has data to construct partial dependence bivariate plot, the first two column includes
#' values of the selected environmental variables, the third column the predicted suitability.
#' \item training_boundaries: has data to plot boundaries of training data.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' # Load data
#' envar <- system.file("external/envar.tif", package = "adm") %>%
#'   rast()
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
#' # Prepare data for Bivariate Partial Dependence Plots
#' bpdp_data <- data_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "sand"),
#'   resolution = 25,
#'   training_data = some_sp,
#'   response_name = "Abundance",
#'   projection_data = envar,
#'   training_boundaries = "convexh"
#' )
#'
#' bpdp_data
#' }
data_abund_bpdp <-
  function(model,
           predictors,
           resolution = 50,
           training_data = NULL,
           invert_transform = NULL,
           response_name = NULL,
           training_boundaries = NULL,
           projection_data = NULL,
           sample_size = NULL,
           training_raster = NULL,
           x_coord = NULL,
           y_coord = NULL) {
    self <- Abundance_inverted <- y <- NULL

    if (!is.null(predictors) & length(predictors) < 2) {
      stop("Please provide at least two predictors.")
    }

    # Extract training data

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

    if (!is.null(training_boundaries) & is.null(training_data)) {
      stop("To extract data to delimit training boundaries it is necessary to provide training data in 'training_data' argument")
    }
    if (!is.null(training_boundaries)) {
      if (!any(training_boundaries %in% c("convexh", "rectangle"))) {
        stop(
          "'training_boundaries' argument could assume one of the following value: NULL, 'convexh', or 'rectangle'"
        )
      }
    }
    if (is.null(training_boundaries)) {
      training_boundaries <- 1
    }

    if (any(class(model)[1] == c("gamlss", "luz_module_fitted", "xgb.Booster"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for GLM, GAM and DNN it is necessary to provide calibration data in 'training_data' argument"
        )
      }
      if (variables[["model"]] == "cnn") {
        x <- training_data[, c(as.vector(variables[1, ])[2:ncol(variables)] %>% unlist(), x_coord, y_coord)]
      } else {
        x <- training_data[, as.vector(variables[1, ])[2:ncol(variables)] %>% unlist()]
      }
    }

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula", "ksvm", "gbm", "maxnet"))) {
      if (is.null(training_data)) {
        stop(
          "To estimate partial plot data for Neural Networks, Random Forest, Support Vector Machine it is necessary to provide calibration data in 'training_data' argument"
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
    if (training_boundaries == "convexh") {
      if (any(sapply(x[predictors], is.factor))) {
        chulld <- NULL
      } else {
        chulld <- x[grDevices::chull(x[predictors]), predictors]
        chulld <- dplyr::as_tibble(chulld)
      }
    } else if (training_boundaries == "rectangle") {
      if (any(sapply(x[predictors], is.factor))) {
        chulld <- NULL
      } else {
        chulld <- apply(x[predictors], 2, range)
        chulld <- expand.grid(chulld[, 1], chulld[, 2], stringsAsFactors = FALSE)
        names(chulld) <- predictors
        chulld <- dplyr::as_tibble(chulld)
      }
    } else {
      chulld <- NULL
    }


    # Control average factor level
    fact <- sapply(x, is.factor)
    suit_c <- which(!fact)
    fact <- which(fact)
    suit_c <- data.frame(t(apply(x[suit_c], 2, mean)))

    if (sum(fact) > 0) {
      for (i in 1:length(fact)) {
        ff <- sort(data.frame(unique(x[names(fact[i])]))[, 1])
        ff <- ff[as.integer(length(ff) / 2)]
        suit_c[names(fact)[i]] <- ff
      }
    }


    if (any(predictors %in% names(fact))) {
      if (is.null(projection_data)) {
        filt <- sapply(x[predictors], is.factor)
        rng1 <- range(x[predictors][!filt])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- sort(data.frame(unique(x[predictors][filt]))[, 1]) # factor
      } else {
        # Range projection data
        projection_data <- projection_data[[predictors]]
        filt <- is.factor(projection_data[[predictors]])
        rng1 <- terra::minmax(projection_data[[!filt]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- as.data.frame(projection_data[[predictors]][[filt]])[, 1] %>% unique()
      }
    } else {
      if (is.null(projection_data)) {
        rng1 <- range(x[, predictors[1]])
        rng2 <- range(x[, predictors[2]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- seq(rng2[1], rng2[2], length.out = resolution)
      } else {
        # Range projection data
        rng1 <- terra::minmax(projection_data[[predictors[1]]])
        rng2 <- terra::minmax(projection_data[[predictors[2]]])
        rng1 <- seq(rng1[1], rng1[2], length.out = resolution)
        rng2 <- seq(rng2[1], rng2[2], length.out = resolution)
      }
    }

    rng <- expand.grid(rng1, rng2, stringsAsFactors = FALSE)
    if (any(sapply(rng, is.factor))) {
      rng <- rng[sapply(rng, is.factor) + 1]
      names(rng) <- predictors[sapply(rng, is.factor) + 1]
    } else {
      names(rng) <- predictors
    }
    suit_c <- suit_c %>% dplyr::select(!{{ predictors }})
    suit_c <- data.frame(rng, suit_c)

    # Make cnn samples
    if (variables[["model"]] == "cnn") {
      suit_c <- suit_c %>%
        dplyr::select(-all_of(x_coord), -all_of(y_coord))

      training_raster <- training_raster[[names(suit_c %>% dplyr::select(-variables[["response"]]))]]

      for (n in names(training_raster)) {
        if (!(n %in% predictors)) {
          training_raster[[n]] <- suit_c[[n]] %>% unique()
        }
      }

      randompoints <- terra::spatSample(training_raster, nrow(rng), xy = T) %>%
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
        random_samples$predictors[[i]][, , predictors[[1]]] <- suit_c[i, predictors][[1]]
        random_samples$predictors[[i]][, , predictors[[2]]] <- suit_c[i, predictors][[2]]
      }
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

      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(
            stats::predict(
              model,
              newdata = pred_dataset,
              type = "response"
            ) %>%
              as.numeric()
          )
        )
    }

    #### dnn ####
    if (class(model)[1] == "luz_module_fitted" & variables[["model"]] == "dnn") {
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

      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(
            stats::predict(
              model,
              newdata = pred_dataset,
              type = "response"
            ) %>%
              as.numeric()
          )
        )
    }

    #### glm and gam ####
    if (class(model)[1] == "gamlss") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(
            predict(
              model,
              newdata = suit_c,
              data = training_data,
              type = "response"
            )
          )
        )
    }

    #### xgb ####
    if (class(model)[1] == "xgb.Booster") {
      pred_matrix <- list(
        data = stats::model.matrix(~ . - 1, data = suit_c[, model$feature_names])
      )

      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(stats::predict(model, newdata = pred_matrix$data, type = "response"))
        )
    }

    #### gbm ####
    if (class(model)[1] == "gbm") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(gbm::predict.gbm(
            model,
            newdata = suit_c, type = "response"
          ))
        )
    }

    #### net ####
    if (class(model)[1] == "nnet.formula") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(
            stats::predict(model, newdata = suit_c, type = "raw")
          )
        )
    }

    #### raf ####
    if (class(model)[1] == "randomForest.formula") {
      suit_c <-
        data.frame(suit_c[1:2], Abundance = suppressMessages(
          stats::predict(model, suit_c, type = "response")
        ))
    }

    #### svm ####
    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = kernlab::predict(model, suit_c, type = "response")
        )
    }

    result <- list("pdpdata" = dplyr::as_tibble(suit_c), "training_boundaries" = chulld)

    if (!is.null(invert_transform)) {
      result[["pdpdata"]] <- result[["pdpdata"]] %>% dplyr::mutate(
        Abundance = adm_transform(result[["pdpdata"]],
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

    if (!is.null(response_name)) {
      idx <- which(names(result[["pdpdata"]]) == "Abundance")
      names(result[["pdpdata"]])[[idx]] <- response_name
    }

    return(result)
  }

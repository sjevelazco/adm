#' Calculate data to construct partial dependence surface plots
#'
#' @param model
#' @param predictors
#' @param resolution
#' @param training_data
#' @param training_boundaries
#' @param projection_data
#'
#' @importFrom dplyr as_tibble select
#' @importFrom gbm predict.gbm
#' @importFrom grDevices chull
#' @importFrom kernlab predict
#' @importFrom stats na.omit
#' @importFrom terra minmax
#' @importFrom torch dataset torch_tensor
#'
#' @return
#' @export
#'
#' @examples
data_abund_bpdp <-
  function(model,
           predictors,
           resolution = 50,
           training_data = NULL,
           invert_transform = NULL,
           response_name = NULL,
           training_boundaries = NULL,
           projection_data = NULL) {
    self <- NULL

    # Extract training data

    if (class(model)[1] == "list") {
      if (all(c("model", "predictors", "performance", "performance_part", "predicted_part") %in% names(model))
      ) {
        variables <- model$predictors
        model <- model[[1]]
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

    if (any(class(model)[1] == c("gamlss", "luz_module_fitted"))) {
      if (is.null(training_data)) {
        stop(
          "For estimating partial plot data for GLM, GAM and DNN it is necessary to provide calibration data in 'training_data' argument"
        )
      }
      x <- training_data[, as.vector(variables[1, ])[2:ncol(variables)] %>% unlist()]
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
        chulld <- expand.grid(chulld[, 1], chulld[, 2])
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

    rng <- expand.grid(rng1, rng2)
    if (any(sapply(rng, is.factor))) {
      rng <- rng[sapply(rng, is.factor) + 1]
      names(rng) <- predictors[sapply(rng, is.factor) + 1]
    } else {
      names(rng) <- predictors
    }
    suit_c <- suit_c %>% dplyr::select(!{{ predictors }})
    suit_c <- data.frame(rng, suit_c)


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

    if (class(model)[1] == "gbm") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(gbm::predict.gbm(
            model,
            newdata = suit_c, type = "response"
          ))
        )
    }


    if (class(model)[1] == "nnet.formula") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = suppressMessages(
            stats::predict(model, newdata = suit_c, type = "raw")
          )
        )
    }

    if (class(model)[1] == "randomForest.formula") {
      suit_c <-
        data.frame(suit_c[1:2], Abundance = suppressMessages(
          stats::predict(model, suit_c, type = "response")
        ))
    }

    if (class(model)[1] == "ksvm") {
      suit_c <-
        data.frame(suit_c[1:2],
          Abundance = kernlab::predict(model, suit_c, type = "response")
        )
    }

    result <- list("pspdata" = dplyr::as_tibble(suit_c), "training_boundaries" = chulld)
    
    if(!is.null(invert_transform)){
      result[["pspdata"]] <- result[["pspdata"]] %>% dplyr::mutate(
        Abundance = adm_transform(result[["pspdata"]],
                                  "Abundance",
                                  invert_transform[["method"]],
                                  inverse = TRUE,
                                  t_terms = c(invert_transform[["a"]] %>% as.numeric,
                                              invert_transform[["b"]] %>% as.numeric)
        ) %>% dplyr::pull(Abundance_inverted)
      )
    }
    
    if(!is.null(response_name)){
          idx <- which(names(result[["pspdata"]])=="Abundance") 
          names(result[["pspdata"]])[[idx]] <- response_name
    }
    
    return(result)
  }

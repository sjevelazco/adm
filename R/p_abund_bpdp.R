#' Bivariate partial dependence plots for abundance-based distribution models
#'
#' @description Create bivariate partial dependence plots to explore the marginal effect of
#' predictors on modeled abundance
#'
#' @param model A model object found in the first element of the list returned
#' by any function from the fit_abund_ or tune_abund_ function families
#' @param predictors character. Vector of predictor name(s) to calculate partial dependence plots.
#' If NULL all predictors will be used. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict abundance
#' values for continuous predictors. Default 50
#' @param training_data data.frame or tibble. Database with response and predictor values used
#' to fit a model. Default NULL
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. When this argument is used, function will calculate partial dependence curves
#' distinguishing conditions used in training and projection conditions
#' (i.e., projection data present in projection area but not training). Default NULL
#' @param training_boundaries character. Plot training conditions boundaries based on training
#' data.
#' If training_boundaries = "convexh", function will delimit training environmental region based on a
#' convex-hull. If training_boundaries = "rectangle", function will delimit training environmental
#' region based on four straight lines. If used any methods it is necessary provide
#' data in training_data argument.If NULL all predictors will be used. Default NULL.
#' @param invert_transform logical. Invert transformation of response variable. Useful for those cases that the response variable was transformed with one of the method in \code{\link{adm_transform}}. Default NULL
#' @param response_name character. Name of the response variable. Default NULL
#' @param color_gradient character. Vector with gradient colors. Default c("#000004", "#1B0A40", "#4A0C69", "#781B6C", "#A42C5F", "#CD4345", "#EC6824", "#FA990B", "#F7CF3D", "#FCFFA4")
#' @param color_training_boundaries character. A vector with one color used to color points of residuals, Default "white"
#' @param set_max numeric. Set a maximum abundance value to plot
#' @param set_min numeric. Set a minimum abundance value to plot
#' @param theme ggplot2 theme. Default ggplot2::theme_classic()
#'
#' @importFrom ggplot2 ggplot aes geom_raster coord_cartesian geom_polygon geom_rect labs scale_fill_gradientn theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stringr str_remove
#' @importFrom utils combn
#'
#' @seealso \code{\link{data_abund_pdp}}, \code{\link{data_abund_bpdp}}, \code{\link{p_abund_pdp}}
#'
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
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
#' # Bivariate Dependence Plots:
#' # In different resolutions
#' p_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "sand"),
#'   training_data = some_sp,
#'   resolution = 50
#' )
#'
#' p_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "sand"),
#'   training_data = some_sp,
#'   resolution = 25
#' )
#'
#' # With projection and training boundaries
#' p_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "elevation", "sand"),
#'   training_data = some_sp,
#'   projection_data = envar,
#'   training_boundaries = "rectangle"
#' )
#'
#' p_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "elevation", "sand"),
#'   training_data = some_sp,
#'   projection_data = envar,
#'   training_boundaries = "convexh"
#' )
#'
#' # Customize colors and theme
#' p_abund_bpdp(
#'   model = mglm,
#'   predictors = c("bio12", "sand"),
#'   training_data = some_sp,
#'   projection_data = envar,
#'   training_boundaries = "convexh",
#'   color_gradient =
#'     c(
#'       "#122414", "#183C26", "#185437", "#106D43", "#0F874C",
#'       "#2D9F54", "#61B463", "#8DC982", "#B3E0A7", "#D7F9D0"
#'     ),
#'   color_training_boundaries = "purple",
#'   theme = ggplot2::theme_dark()
#' )
#' }
p_abund_bpdp <-
  function(model,
           predictors = NULL,
           resolution = 50,
           training_data = NULL,
           projection_data = NULL,
           training_boundaries = NULL,
           invert_transform = NULL,
           response_name = NULL,
           color_gradient = c(
             "#000004",
             "#1B0A40",
             "#4A0C69",
             "#781B6C",
             "#A42C5F",
             "#CD4345",
             "#EC6824",
             "#FA990B",
             "#F7CF3D",
             "#FCFFA4"
           ),
           color_training_boundaries = "white",
           set_max = NULL,
           set_min = NULL,
           theme = ggplot2::theme_classic()) {
    Abundance <- Type <- Value <- val <- sym <- NULL

    if (!is.null(predictors) & length(predictors) < 2) {
      stop("Please provide at least two predictors.")
    }

    if (class(model)[1] == "list") {
      if (all(c("model", "predictors", "performance", "performance_part", "predicted_part") %in% names(model))
      ) {
        variables <- model$predictors
        model_l <- model
        model <- model[[1]]
      }
    } else {
      stop('Please, use tune_abund_ or fit_abund_ output list in "model" argument.')
    }

    if (!is.null(training_boundaries) & is.null(training_data)) {
      stop(
        "To plot bivariate partial dependence plot with training condition boundaries it is necessary to provide calibration data in 'training_data' argument"
      )
    }
    if (!is.null(training_boundaries)) {
      if (!any(training_boundaries %in% c("convexh", "rectangle"))) {
        stop(
          "'training_boundaries' argument could assume one of the following value: NULL, 'convexh', or 'rectangle'"
        )
      }
    }

    if (class(model)[1] %in% c("luz_module_fitted", "xgb.Booster")) {
      if (!is.null(training_data)) {
        v <- training_data[variables[1, 3:ncol(variables)] %>%
          as.vector() %>%
          unlist()] %>% sapply(class)
      } else {
        stop("Training data needed.")
      }
    }

    if (class(model)[1] == "gamlss") {
      if (variables[["model"]] == "gam") {
        v <- attr(model$mu.terms, "dataClasses")[-1]
        names(v) <- names(v) %>%
          stringr::str_remove("pb\\(") %>%
          stringr::str_remove("\\)")
      } else if (variables[["model"]] == "glm") {
        v <- attr(model$mu.terms, "dataClasses")[-1]
      }
    }


    if (class(model)[1] == "gbm") {
      v <- attr(model$Terms, "dataClasses")[-1]
    }

    if (any(class(model)[1] == c("nnet.formula", "randomForest.formula"))) {
      v <- attr(model$terms, "dataClasses")[-1]
    }

    if (class(model)[1] == "ksvm") {
      v <- attr(model@terms, "dataClasses")[-1]
    }

    v <- v[order(names(v))]
    if (!is.null(predictors)) {
      v <- v[names(v) %in% predictors]
    }

    var_comb <- utils::combn(names(v), 2) %>%
      t() %>%
      data.frame()

    if (any(v == "factor")) {
      filt <- (v[var_comb[, 1]] == "factor" & v[var_comb[, 2]] == "factor")
      var_comb <- var_comb[!filt, ]
      filt <- which(v[var_comb[, 1]] == "factor")
      if (length(filt) > 0) {
        var_comb[filt, ] <- var_comb[filt, 2:1]
      }
    }


    if (is.null(response_name)) {
      response_name <- "Abundance"
    }

    p_list <- list()
    abundance_values <- c()

    for (i in 1:nrow(var_comb)) {
      xenv <- var_comb[i, 1]
      yenv <- var_comb[i, 2]

      crv <-
        data_abund_bpdp(
          model = model_l,
          predictors = c(xenv, yenv),
          resolution = resolution,
          training_boundaries = training_boundaries,
          projection_data = projection_data,
          training_data = training_data,
          response_name = response_name,
          invert_transform = invert_transform
        )

      # Coleta os valores de abund<U+00E2>ncia para calcular min e max
      if (!is.null(set_max)) {
        if (is.numeric(set_max)) {
          crv[[1]][which(crv[[1]] %>% dplyr::pull(response_name) > set_max), response_name] <- set_max
        } else {
          stop("set_max and set_min should be numerical.")
        }
      }

      if (!is.null(set_min)) {
        if (is.numeric(set_min)) {
          crv[[1]][which(crv[[1]] %>% dplyr::pull(response_name) < set_min), response_name] <- set_min
        } else {
          stop("set_max and set_min should be numerical.")
        }
      }

      abundance_values <- c(abundance_values, crv[[1]][[response_name]])

      v1 <- names(crv[[1]])[1]
      v2 <- names(crv[[1]])[2]
      names(crv[[1]])[1:2] <- c("v1", "v2")

      p_list[[i]] <-
        ggplot2::ggplot(crv[[1]], ggplot2::aes(v1, v2)) +
        ggplot2::geom_raster(ggplot2::aes(fill = !!sym(response_name))) +
        ggplot2::coord_cartesian(expand = FALSE) +
        {
          if (!is.null(training_boundaries) & !is.null(crv$training_boundaries)) {
            if (training_boundaries == "convexh") {
              names(crv[[2]]) <- c("v1", "v2")
              ggplot2::geom_polygon(
                data = crv[[2]],
                ggplot2::aes(v1, v2),
                color = color_training_boundaries,
                fill = "transparent"
              )
            }
          }
        } +
        {
          if (!is.null(training_boundaries) & !is.null(crv$training_boundaries)) {
            if (training_boundaries == "rectangle") {
              names(crv[[2]]) <- c("v1", "v2")
              ggplot2::geom_rect(
                data = crv[[2]],
                ggplot2::aes(
                  xmin = min(v1), xmax = max(v1),
                  ymin = min(v2), ymax = max(v2)
                ),
                color = color_training_boundaries,
                fill = "transparent"
              )
            }
          }
        } +
        ggplot2::labs(x = v1, y = v2)
    }


    max_abund <- max(abundance_values, na.rm = TRUE)
    min_abund <- min(abundance_values, na.rm = TRUE)


    for (i in 1:length(p_list)) {
      p_list[[i]] <- p_list[[i]] +
        ggplot2::scale_fill_gradientn(colours = color_gradient, limits = c(min_abund, max_abund)) +
        theme
    }

    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      result <- patchwork::wrap_plots(p_list) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::theme(legend.position = "bottom") +
          theme(legend.title = element_blank())
      return(result)
    }
  }

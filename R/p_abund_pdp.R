#' Partial dependent plots for abundance-based distribution models
#'
#' @description Create partial dependence plots to explore the marginal effect of
#' predictors on modeled abundance
#'
#' @param model A model object found in the first element of the list returned
#' by any function from the fit_abund_ or tune_abund_ function families
#' @param predictors character. Vector of predictor name(s) to calculate partial dependence plots.
#' If NULL all predictors will be used. Default NULL
#' @param resolution numeric. Number of equally spaced points at which to predict abundance
#' values for continuous predictors. Default 50
#' @param resid logical. Calculate residuals based on training data. Default FALSE
#' @param training_data data.frame or tibble. Database with response and predictor values used
#' to fit a model. Default NULL
#' @param projection_data SpatRaster. Raster layer with environmental variables used for model
#' projection. When this argument is used, function will calculate partial dependence curves
#' distinguishing conditions used in training and projection conditions
#' (i.e., projection data present in projection area but not training). Default NULL
#' @param rug logical. Add rug plot to partial dependence plot. Default FALSE
#' @param colorl character. Vector with colors to plot partial dependence curves. Default c("#462777", "#6DCC57")
#' @param colorp character. Color to plot residuals. Default "black"
#' @param alpha numeric. Transparency of residuals. Default 0.2
#' @param theme ggplot2 theme. Default ggplot2::theme_classic()
#' @param invert_transform logical. If TRUE, inverse transformation of response variable will be applied.
#' @param response_name character. Name of the response variable. Default NULL
#'
#' @details This function creates partial dependent plots to explore the marginal effect of
#' predictors on modeled abundance. If projection_data is used, function will extract the minimum and
#' maximum values found in a region or time period to which a model will be projected.
#' If the range of projection data is greater than of the training data it will be
#' plotted with a different color. Partial dependence plot could be used to interpret a
#' model or to explore how a model may extrapolate outside the environmental conditions
#' used to train the model.
#'
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot aes scale_y_continuous labs geom_point geom_line geom_rug geom_col scale_color_manual geom_vline theme element_blank
#' @importFrom patchwork wrap_plots plot_layout
#'
#' @seealso \code{\link{data_abund_pdp}}, \code{\link{data_abund_bpdp}}, \code{\link{p_abund_bpdp}}
#'
#' @return A ggplot object
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
#' # Partial Dependence Plots:
#'
#' # In different resolutions
#' p_abund_pdp(
#'   model = mglm,
#'   resolution = 50,
#'   training_data = some_sp,
#'   response_name = "Abundance"
#' )
#'
#' p_abund_pdp(
#'   model = mglm,
#'   resolution = 5,
#'   training_data = some_sp,
#'   response_name = "Abundance"
#' )
#'
#' # Especific variables and different resulotions
#' p_abund_pdp(
#'   model = mglm,
#'   predictors = c("bio12", "sand"),
#'   training_data = some_sp,
#'   response_name = "Abundance"
#' )
#'
#' # With residuals and rug plot
#' p_abund_pdp(
#'   model = mglm,
#'   training_data = some_sp,
#'   response_name = "Abundance",
#'   resid = TRUE
#' )
#'
#' p_abund_pdp(
#'   model = mglm,
#'   training_data = some_sp,
#'   response_name = "Abundance",
#'   rug = TRUE
#' )
#'
#' p_abund_pdp(
#'   model = mglm,
#'   training_data = some_sp,
#'   response_name = "Abundance",
#'   resid = TRUE,
#'   rug = TRUE
#' )
#'
#' # Partial depence plot for training and projection condition found in a projection area
#' p_abund_pdp(
#'   model = mglm,
#'   training_data = some_sp,
#'   projection_data = envar,
#'   response_name = "Abundance",
#'   rug = TRUE
#' )
#'
#' # Custumize colors and theme
#' p_abund_pdp(
#'   model = mglm,
#'   predictors = NULL,
#'   resolution = 100,
#'   resid = TRUE,
#'   training_data = some_sp,
#'   projection_data = envar,
#'   colorl = c("blue", "red"),
#'   colorp = "darkgray",
#'   alpha = 0.4,
#'   theme = ggplot2::theme_dark()
#' )
#' }
p_abund_pdp <-
  function(model,
           predictors = NULL,
           resolution = 100,
           resid = FALSE,
           training_data = NULL,
           invert_transform = NULL,
           response_name = NULL,
           projection_data = NULL,
           rug = FALSE,
           colorl = c("#462777", "#6DCC57"),
           colorp = "black",
           alpha = 0.2,
           theme = ggplot2::theme_classic()) {
    Type <- Value <- val <- Abundance <- sym <- NULL

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

    if (!all(variables[1, 2:ncol(variables)] %>%
      as.vector() %>%
      unlist() %in% names(training_data))) {
      stop("Variables not present in training data. Did you use the wrong dataset?")
    }

    if (class(model)[1] %in% c("luz_module_fitted", "xgb.Booster")) {
      if (!is.null(training_data)) {
        v <- training_data[variables[1, 2:ncol(variables)] %>%
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
        v <- v[names(v) %in% variables[1, 2:ncol(variables)] %>%
          as.vector() %>%
          unlist()]
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

    if (is.null(response_name)) {
      response_name <- "Abundance"
    }

    p <- list()

    if (is.null(projection_data)) {
      for (i in 1:length(v)) {
        crv <-
          data_abund_pdp(
            model = model_l,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = NULL,
            training_data = training_data,
            invert_transform = invert_transform,
            response_name = response_name
          )

        if (v[i] == "numeric") {
          xn <- data.frame(crv[[1]])[, 1]
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(x = !!xn, y = !!sym(response_name))) +
            ggplot2::labs(x = names(crv[[1]])[1]) +
            {
              if (resid) {
                xn2 <- data.frame(crv[[2]])[, 1]
                ggplot2::geom_point(
                  data = crv[[2]], color = colorp,
                  ggplot2::aes(!!xn2, !!sym(response_name)), alpha = alpha
                )
              }
            } +
            ggplot2::geom_line(col = rev(colorl)[1], linewidth = 0.8)

          if (rug) {
            xn2 <- data.frame(crv[[2]])[, 1]

            rug_df <- crv[[2]] %>% mutate(!!response_name := max(crv[[1]][response_name]))
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = rug_df,
                ggplot2::aes(!!xn2, !!sym(response_name)),
                sides = "b",
                alpha = 0.3
              )
          }
        } else {
          xn <- data.frame(crv[[1]])[, 1]
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, !!sym(response_name))) +
            ggplot2::geom_col(fill = rev(colorl)[1]) +
            ggplot2::labs(x = names(crv[[1]])[1])
        }
      }
    } else {
      for (i in 1:length(v)) {
        crv <-
          data_abund_pdp(
            model = model_l,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = projection_data[[c(names(v[i]))]],
            training_data = training_data,
            invert_transform = invert_transform,
            response_name = response_name
          )

        if (v[i] == "numeric") {
          rvar <- range(crv[[1]][crv[[1]]$Type == "Training", names(v[i])])
          xn <- data.frame(crv[[1]])[, 1]

          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, !!sym(response_name))) +
            ggplot2::labs(x = names(crv[[1]])[1]) +
            {
              if (resid) {
                xn2 <- data.frame(crv[[2]])[, 1]
                ggplot2::geom_point(
                  data = crv[[2]], ggplot2::aes(!!xn2, !!sym(response_name)),
                  alpha = alpha, color = colorp
                )
              }
            } +
            ggplot2::geom_line(ggplot2::aes(color = Type, group = 1), linewidth = 0.8) +
            ggplot2::scale_color_manual(
              values = colorl,
              breaks = c("Projection", "Training"),
              name = "Range"
            ) +
            ggplot2::geom_vline(
              xintercept = rvar,
              col = "gray70",
              linetype = 2
            )

          if (rug) {
            xn2 <- data.frame(crv[[2]])[, 1]

            rug_df <- crv[[2]] %>% mutate(!!response_name := max(crv[[1]][response_name]))
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = rug_df,
                ggplot2::aes(!!xn2, !!sym(response_name)),
                sides = "b",
                alpha = 0.5
              )
          }
        } else {
          xn <- crv[[1]] %>% dplyr::pull(names(crv[[1]])[1])
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, !!sym(response_name))) +
            ggplot2::geom_col(fill = rev(colorl)[1]) +
            ggplot2::labs(x = names(crv[[1]])[1])
        }
      }
    }

    # Theme
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme
    }

    # Remove y axis titles
    if (length(p) >= 4) {
      sq <- length(p) / round(sqrt(length(p)))
      sq <- seq(1, length(p), by = sq)
      sq2 <- 1:length(p)
      sq2 <- sq2[!sq2 %in% sq]
    } else if (length(p) < 4 & length(p) > 2) {
      sq2 <- 2:length(p)
    }
    if (exists("sq2")) {
      for (i in sq2) {
        p[[i]] <- p[[i]] + ggplot2::theme(axis.title.y = ggplot2::element_blank())
      }
    }

    # ncol = round(sqrt(length(p)))
    patchwork::wrap_plots(p) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank())
  }

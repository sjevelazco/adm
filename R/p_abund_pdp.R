#' ADM Partial Dependent Plot
#'
#' @param model
#' @param predictors
#' @param resolution
#' @param resid
#' @param training_data
#' @param projection_data
#' @param clamping
#' @param rug
#' @param colorl
#' @param colorp
#' @param alpha
#' @param theme
#'
#' @importFrom dplyr pull
#' @importFrom ggplot2 ggplot aes scale_y_continuous labs geom_point geom_line geom_rug geom_col scale_color_manual geom_vline theme element_blank
#' @importFrom patchwork wrap_plots plot_layout
#'
#' @return
#' @export
#'
#' @examples
p_abund_pdp <-
  function(model,
           predictors = NULL,
           resolution = 100,
           resid = FALSE,
           training_data = NULL,
           projection_data = NULL,
           clamping = FALSE,
           rug = FALSE,
           colorl = c("#462777", "#6DCC57"),
           colorp = "black",
           alpha = 0.2,
           theme = ggplot2::theme_classic()) {
    Type <- Value <- val <- Abundance <- NULL

    # TODO
    # if (class(model)[1] == "gam") {
    #   v <- attr(model$terms, "dataClasses")[-1]
    # }

    # TODO
    # if (class(model)[1] == "graf") {
    #   v <- sapply(model$obsx, class)
    # }

    # TODO
    # if (class(model)[1] == "glm") {
    #   flt <- grepl("[I(]", attr(model$terms, "term.labels")) |
    #     grepl(":", attr(model$terms, "term.labels"))
    #   flt <- attr(model$terms, "term.labels")[!flt]
    #   v <- attr(model$terms, "dataClasses")[flt]
    # }

    # TODO
    # if (class(model)[1] == "xgb.Booster"){
    #   if (!is.null(training_data)){
    #     v <- training_data[model$feature_names] %>% sapply(class)
    #   } else {
    #     stop("Training data needed.")
    #   }
    # }

    if (class(model)[1] == "list") {
      if (all(c("model", "predictors", "performance", "performance_part", "predicted_part") %in% names(model))
      ) {
        variables <- model$predictors
        model <- model[[1]]
      }
    }

    if(!all(variables[1, 2:ncol(variables)] %>%
            as.vector() %>%
            unlist() %in% names(training_data))){
      stop("Variables not present in training data. Did you use the wrong dataset?")
    }
    
    if (class(model)[1] == "luz_module_fitted") {
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
        v <- v[names(v)%in%variables[1, 2:ncol(variables)] %>%
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

    p <- list()

    if (is.null(projection_data)) {
      for (i in 1:length(v)) {
        crv <-
          data_abund_pdp(
            model = model,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = NULL,
            training_data = training_data
          )

        if (v[i] == "numeric") {
          xn <- data.frame(crv[[1]])[, 1]
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(x = !!xn, y = Abundance)) +
            # ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::labs(x = names(crv[[1]])[1]) +
            {
              if (resid) {
                xn2 <- data.frame(crv[[2]])[, 1]
                ggplot2::geom_point(
                  data = crv[[2]], color = colorp,
                  ggplot2::aes(!!xn2, Abundance), alpha = alpha
                )
              }
            } +
            ggplot2::geom_line(col = rev(colorl)[1], linewidth = 0.8)

          if (rug) {
            xn2 <- data.frame(crv[[2]])[, 1]
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = crv[[2]],
                ggplot2::aes(!!xn2, Abundance),
                sides = "b",
                alpha = 0.3
              )
          }
        } else {
          xn <- data.frame(crv[[1]])[, 1]
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, Abundance)) +
            # ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::geom_col(fill = rev(colorl)[1]) +
            ggplot2::labs(x = names(crv[[1]])[1])
        }
      }
    } else {
      for (i in 1:length(v)) {
        crv <-
          data_abund_pdp(
            model = model,
            predictors = names(v[i]),
            resolution = resolution,
            resid = any(c(resid, rug)),
            projection_data = projection_data[[c(names(v[i]))]],
            training_data = training_data
          )

        if (v[i] == "numeric") {
          rvar <- range(crv[[1]][crv[[1]]$Type == "Training", names(v[i])])
          xn <- data.frame(crv[[1]])[, 1]

          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, Abundance)) +
            ggplot2::labs(x = names(crv[[1]])[1]) +
            {
              if (resid) {
                xn2 <- data.frame(crv[[2]])[, 1]
                ggplot2::geom_point(
                  data = crv[[2]], ggplot2::aes(!!xn2, Abundance),
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
            # ggplot2::scale_y_continuous(limits = c(0, 1)) +
            ggplot2::geom_vline(
              xintercept = rvar,
              col = "gray70",
              linetype = 2
            )

          if (rug) {
            xn2 <- data.frame(crv[[2]])[, 1]
            p[[i]] <- p[[i]] +
              ggplot2::geom_rug(
                data = crv[[2]],
                ggplot2::aes(!!xn2, Abundance),
                sides = "b",
                alpha = 0.5
              )
          }
        } else {
          xn <- crv[[1]] %>% dplyr::pull(names(crv[[1]])[1])
          p[[i]] <-
            ggplot2::ggplot(crv[[1]], ggplot2::aes(!!xn, Abundance)) +
            # ggplot2::scale_y_continuous(limits = c(0, 1)) +
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

#' ADM Bivariate partial dependence plot
#'
#' @param model 
#' @param predictors 
#' @param resolution 
#' @param training_data 
#' @param training_boundaries 
#' @param projection_data 
#' @param clamping 
#' @param color_gradient 
#' @param color_training_boundaries 
#' @param theme 
#' 
#' @importFrom ggplot2 ggplot aes geom_raster coord_cartesian geom_polygon geom_rect labs scale_fill_gradientn theme
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stringr str_remove
#' @importFrom utils combn
#'
#' @return
#' @export
#'
#' @examples
p_abund_bpdp <-
  function(model,
           predictors = NULL,
           resolution = 50,
           training_data = NULL,
           training_boundaries = NULL,
           projection_data = NULL,
           clamping = FALSE,
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
           theme = ggplot2::theme_classic()) {
    Abundance <- Type <- Value <- val <- NULL
    
    if (class(model)[1]=="list"){
      if(all(names(model) %in% c("model","predictors","performance","performance_part","predicted_part"))
      ){      
        variables <- model$predictors 
        model <- model[[1]]
      }
    }
    
    if (!is.null(training_boundaries) & is.null(training_data)) {
      stop(
        "To plot bivariate partial dependence plot with training condition boundaries it is necessary to provide calibration data in 'training_data' argument"
      )
    }
    if(!is.null(training_boundaries)){
      if(!any(training_boundaries %in% c("convexh", "rectangle"))){
        stop(
          "'training_boundaries' argument could assume one of the following value: NULL, 'convexh', or 'rectangle'"
        )
      }
    }
    
    if (class(model)[1] == "luz_module_fitted"){
      if (!is.null(training_data)){
        v <- training_data[variables[1,2:ncol(variables)]%>%as.vector%>%unlist] %>% sapply(class)
      } else {
        stop("Training data needed.")
      }
    }
    
    if (class(model)[1] == "gamlss") {
      if (variables[["model"]] == "gam"){
        v <- attr(model$mu.terms, "dataClasses")[-1]
        names(v) <- names(v) %>% stringr::str_remove("pb\\(") %>% stringr::str_remove("\\)")
      } else if (variables[["model"]] == "glm"){
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
    
    p_list <- list()
    abundance_values <- c()
    
    for (i in 1:nrow(var_comb)) {
      xenv <- var_comb[i, 1]
      yenv <- var_comb[i, 2]
      
      crv <-
        data_abund_bpdp(
          model = model,
          predictors = c(xenv, yenv),
          resolution = resolution,
          training_boundaries = training_boundaries,
          projection_data = projection_data,
          training_data = training_data,
          clamping = clamping
        )
      
      # Coleta os valores de abundância para calcular min e max
      abundance_values <- c(abundance_values, crv[[1]]$Abundance)
      
      v1 <- names(crv[[1]])[1]
      v2 <- names(crv[[1]])[2]
      names(crv[[1]])[1:2] <- c("v1", "v2")
      
      p_list[[i]] <-
        ggplot2::ggplot(crv[[1]], ggplot2::aes(v1, v2)) +
        ggplot2::geom_raster(aes(fill = Abundance)) +
        ggplot2::coord_cartesian(expand = FALSE) +
        {
          if(!is.null(training_boundaries) & !is.null(crv$training_boundaries)){
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
          if(!is.null(training_boundaries) & !is.null(crv$training_boundaries)){
            if (training_boundaries == "rectangle") {
              names(crv[[2]]) <- c("v1", "v2")
              ggplot2::geom_rect(
                data = crv[[2]],
                ggplot2::aes(xmin = min(v1), xmax = max(v1),
                             ymin = min(v2), ymax = max(v2)),
                color = color_training_boundaries,
                fill = "transparent"
              )
            }
          }
        } +
        ggplot2::labs(x = v1, y = v2)
    }
    
   
    min_abund <- min(abundance_values, na.rm = TRUE)
    max_abund <- max(abundance_values, na.rm = TRUE)
    
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

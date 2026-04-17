#' Calculate abundance distribution model uncertainty
#'
#' @description This function calculates the uncertainty of an abundance distribution model by performing a bootstrap procedure.
#' It refits the model multiple times on resampled data and then calculates the standard deviation of the predictions across
#' all iterations.
#'
#' @param models A model object from `fit_abund_*` or `tune_abund_*` functions.
#' @param training_data A data.frame or tibble with abundance data and predictors.
#' @param response character. Column name of the response variable.
#' @param projection_data A SpatRaster object with the environmental layers for projection.
#' @param iteration numeric. The number of bootstrap iterations. Default 50.
#' @param n_cores numeric. The number of cores to use for parallel processing. Default 1.
#' @param ... Additional arguments passed to refitting functions or \code{\link{adm_predict}} 
#' (e.g., \code{x}, \code{y}, \code{rasters}, \code{sample_size} for CNN; \code{custom_architecture} for DNN/CNN; 
#' \code{invert_transform}, \code{transform_negative} for spatial prediction).
#'
#' @return A SpatRaster object with a single layer representing the model uncertainty, calculated as the standard deviation of the bootstrap predictions.
#' @export
#'
#' @seealso \code{\link{adm_predict}}
#'
#' @examples
#' \dontrun{
#' require(terra)
#' require(dplyr)
#'
#' # Load data
#' data("cretusa_data")
#' cretusa_predictors <- system.file("external/cretusa_predictors.tif", package = "adm")
#' cretusa_predictors <- terra::rast(cretusa_predictors)
#'
#' species_data <- adm_extract(
#'   data = cretusa_data,
#'   x = "x",
#'   y = "y",
#'   env_layer = cretusa_predictors
#' )
#'
#' # Fit model
#' mraf <- fit_abund_raf(
#'   data = species_data,
#'   response = "ind_ha",
#'   predictors = c("PC1", "PC2", "PC3"),
#'   partition = ".part"
#' )
#'
#' # Calculate uncertainty
#' unc <- adm_uncertainty(
#'   models = mraf,
#'   training_data = species_data,
#'   response = "ind_ha",
#'   projection_data = cretusa_predictors[[c("PC1", "PC2", "PC3")]],
#'   iteration = 10,
#'   n_cores = 2
#' )
#'
#' plot(unc)
#' }
adm_uncertainty <- function(
  models,
  training_data,
  response,
  projection_data,
  iteration = 50,
  n_cores = 1,
  ...
) {
  # Extract algorithm type
  clss <- models$predictors$model

  # Prepare training dataset (remove partition columns)
  training_data_clean <- training_data %>% dplyr::select(-dplyr::starts_with(".part"))

  # Predictor names
  pr_c <- models$predictors %>%
    dplyr::select(dplyr::starts_with("c")) %>%
    unlist() %>% as.character()
  pr_f <- models$predictors %>%
    dplyr::select(dplyr::starts_with("f")) %>%
    unlist() %>% as.character()
  names(pr_c) <- NULL
  names(pr_f) <- NULL

  # Capture extra arguments for refitting and prediction
  extra_args <- list(...)

  #### Bootstrap approach ####
  my_cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(my_cluster)

  # Refitting and prediction loop
  r_list <- foreach::foreach(
    ii = 1:iteration,
    .packages = c("adm", "dplyr", "terra"),
    .errorhandling = "pass"
  ) %dopar% {
    set.seed(ii)
    
    # Bootstrap sample
    db <- training_data_clean[sample(nrow(training_data_clean), replace = TRUE), ]

    # Refit model based on algorithm type
    m_refit <- switch(clss,
      "raf" = {
        fit_abund_raf(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      mtry = if(!is.null(models$optimal_combination$mtry)) models$optimal_combination$mtry else models$model$mtry,
                      ntree = if(!is.null(models$optimal_combination$ntree)) models$optimal_combination$ntree else models$model$ntree,
                      partition = NULL, verbose = FALSE)
      },
      "glm" = {
        fit_abund_glm(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      distribution = if(!is.null(models$optimal_combination$distribution)) models$optimal_combination$distribution else models$model$family[1],
                      poly = if(!is.null(models$optimal_combination$poly)) models$optimal_combination$poly else 0,
                      inter_order = if(!is.null(models$optimal_combination$inter_order)) models$optimal_combination$inter_order else 0,
                      partition = NULL, verbose = FALSE)
      },
      "gam" = {
        fit_abund_gam(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      distribution = if(!is.null(models$optimal_combination$distribution)) models$optimal_combination$distribution else models$model$family[1],
                      inter = if(!is.null(models$optimal_combination$inter)) models$optimal_combination$inter else "automatic",
                      partition = NULL, verbose = FALSE)
      },
      "gbm" = {
        fit_abund_gbm(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      distribution = if(!is.null(models$optimal_combination$distribution)) models$optimal_combination$distribution else models$model$distribution,
                      n.trees = if(!is.null(models$optimal_combination$n.trees)) models$optimal_combination$n.trees else models$model$n.trees,
                      interaction.depth = if(!is.null(models$optimal_combination$interaction.depth)) models$optimal_combination$interaction.depth else models$model$interaction.depth,
                      n.minobsinnode = if(!is.null(models$optimal_combination$n.minobsinnode)) models$optimal_combination$n.minobsinnode else models$model$n.minobsinnode,
                      shrinkage = if(!is.null(models$optimal_combination$shrinkage)) models$optimal_combination$shrinkage else models$model$shrinkage,
                      partition = NULL, verbose = FALSE)
      },
      "svm" = {
        fit_abund_svm(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      kernel = if(!is.null(models$optimal_combination$kernel)) models$optimal_combination$kernel else "rbfdot",
                      sigma = if(!is.null(models$optimal_combination$sigma)) models$optimal_combination$sigma else "automatic",
                      C = if(!is.null(models$optimal_combination$C)) models$optimal_combination$C else 1,
                      partition = NULL, verbose = FALSE)
      },
      "xgb" = {
        fit_abund_xgb(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      nrounds = if(!is.null(models$optimal_combination$nrounds)) models$optimal_combination$nrounds else 100,
                      max_depth = if(!is.null(models$optimal_combination$max_depth)) models$optimal_combination$max_depth else 5,
                      eta = if(!is.null(models$optimal_combination$eta)) models$optimal_combination$eta else 0.1,
                      gamma = if(!is.null(models$optimal_combination$gamma)) models$optimal_combination$gamma else 1,
                      colsample_bytree = if(!is.null(models$optimal_combination$colsample_bytree)) models$optimal_combination$colsample_bytree else 1,
                      min_child_weight = if(!is.null(models$optimal_combination$min_child_weight)) models$optimal_combination$min_child_weight else 1,
                      subsample = if(!is.null(models$optimal_combination$subsample)) models$optimal_combination$subsample else 0.5,
                      objective = if(!is.null(models$optimal_combination$objective)) models$optimal_combination$objective else "reg:squarederror",
                      partition = NULL, verbose = FALSE)
      },
      "net" = {
        fit_abund_net(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      size = if(!is.null(models$optimal_combination$size)) models$optimal_combination$size else 10,
                      decay = if(!is.null(models$optimal_combination$decay)) models$optimal_combination$decay else 0.1,
                      partition = NULL, verbose = FALSE)
      },
      "dnn" = {
        fit_abund_dnn(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      learning_rate = if(!is.null(models$optimal_combination$learning_rate)) models$optimal_combination$learning_rate else 0.01,
                      n_epochs = if(!is.null(models$optimal_combination$n_epochs)) models$optimal_combination$n_epochs else 10,
                      batch_size = if(!is.null(models$optimal_combination$batch_size)) models$optimal_combination$batch_size else 32,
                      custom_architecture = extra_args$custom_architecture,
                      partition = NULL, verbose = FALSE)
      },
      "cnn" = {
        fit_abund_cnn(data = db, response = response, predictors = pr_c, predictors_f = pr_f,
                      x = extra_args$x, y = extra_args$y, rasters = extra_args$rasters, 
                      sample_size = extra_args$sample_size,
                      learning_rate = if(!is.null(models$optimal_combination$learning_rate)) models$optimal_combination$learning_rate else 0.01,
                      n_epochs = if(!is.null(models$optimal_combination$n_epochs)) models$optimal_combination$n_epochs else 10,
                      batch_size = if(!is.null(models$optimal_combination$batch_size)) models$optimal_combination$batch_size else 32,
                      custom_architecture = extra_args$custom_architecture,
                      partition = NULL, verbose = FALSE)
      }
    )

    # Predict
    p <- adm_predict(models = m_refit, pred = projection_data, training_data = db, ...)
    
    # Return cell values as vector
    as.vector(p[[1]])
  }
  parallel::stopCluster(my_cluster)

  # Identify successful iterations
  r_list_clean <- Filter(function(x) is.numeric(x), r_list)
  if (length(r_list_clean) == 0) stop("All bootstrap iterations failed.")
  
  # Calculate cell-wise standard deviation
  res_mat <- do.call(cbind, r_list_clean)
  unc_vals <- apply(res_mat, 1, stats::sd, na.rm = TRUE)
  
  # Map back to SpatRaster
  res_raster <- projection_data[[1]]
  res_raster[] <- unc_vals
  
  names(res_raster) <- "uncertainty"
  return(res_raster)
}

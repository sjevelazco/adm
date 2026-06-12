#' Calculate abundance distribution model uncertainty
#'
#' @description This function calculates the uncertainty of an abundance distribution model by performing a bootstrap procedure.
#' It refits the model multiple times on resampled data and then calculates the standard deviation of the predictions across
#' all iterations.
#'
#' @param models A model object from `fit_abund_*` or `tune_abund_*` functions.
#' @param training_data A data.frame or tibble with abundance data and predictors.
#' @param response character. Column name of the response variable.
#' @param pred A SpatRaster object with the environmental layers for projection.
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
#'   pred = cretusa_predictors[[c("PC1", "PC2", "PC3")]],
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
  pred,
  iteration = 50,
  sample_prop = 0.80,
  n_cores = 1
) {
  # Extract algorithm type
  clss <- models$predictors$model

  # Predictor names
  pr_c <- models$predictors %>%
    dplyr::select(dplyr::starts_with("c")) %>%
    unlist()
  pr_f <- models$predictors %>%
    dplyr::select(dplyr::starts_with("f")) %>%
    unlist()
  names(pr_c) <- NULL
  names(pr_f) <- NULL

  # #### Bootstrap approach ####
  # my_cluster <- parallel::makeCluster(n_cores)
  # doParallel::registerDoParallel(my_cluster)

  # # Refitting and prediction loop
  # r_list <- foreach::foreach(
  #   ii = 1:iteration,
  #   .packages = c("dplyr", "terra"),
  #   .export = c(
  #     "fit_abund_cnn", "fit_abund_dnn", "fit_abund_gam", "fit_abund_glm",
  #     "fit_abund_net", "fit_abund_qrf", "fit_abund_raf", "fit_abund_svm",
  #     "fit_abund_xgb", "adm_predict"
  #   ),
  #   .errorhandling = "pass"
  # ) %dopar% {
  r_list <- list()
  for (ii in 1:iteration) {
    set.seed(ii)

    # Bootstrap sample
    db <- training_data %>%
      mutate(
        pr_ab = case_when(
          .data[[response]] == 0 ~ 0,
          .data[[response]] > 0 ~ 1
        )
      ) %>%
      group_by(pr_ab) %>%
      slice_sample(prop = sample_prop) %>%
      ungroup()

    # Family for GAM and GLM
    if (clss == "gam" || clss == "glm") {
      fam_char <- if (!is.null(models$optimal_combination$distribution)) models$optimal_combination$distribution else models$model$family[[1]]
      fam_char <- as.character(fam_char)
    }


    # Refit model based on algorithm type
    m_refit <- switch(clss,
      "raf" = {
        fit_abund_raf(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          mtry = models$metadata$hyperparameters$mtry,
          ntree = models$metadata$hyperparameters$ntree,
          partition = NULL, verbose = FALSE
        )
      },
      "glm" = {
        fit_abund_glm(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          distribution = models$metadata$hyperparameters$distribution,
          poly = models$metadata$hyperparameters$poly,
          inter_order = models$metadata$hyperparameters$inter_order,
          partition = NULL, verbose = FALSE
        )
      },
      "gam" = {
        fit_abund_gam(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          inter = models$metadata$hyperparameters$inter,
          distribution = models$metadata$hyperparameters$distribution,
          partition = NULL, verbose = FALSE
        )
      },
      "gbm" = {
        fit_abund_gbm(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          distribution = models$metadata$hyperparameters$distribution,
          n.trees = models$metadata$hyperparameters$n.trees,
          interaction.depth = models$metadata$hyperparameters$interaction.depth,
          n.minobsinnode = models$metadata$hyperparameters$n.minobsinnode,
          shrinkage = models$metadata$hyperparameters$shrinkage,
          partition = NULL, verbose = FALSE
        )
      },
      "svm" = {
        fit_abund_svm(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          kernel = models$metadata$hyperparameters$kernel,
          sigma = models$metadata$hyperparameters$sigma,
          C = models$metadata$hyperparameters$C,
          partition = NULL, verbose = FALSE
        )
      },
      "xgb" = {
        fit_abund_xgb(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          nrounds = models$metadata$boosted_rounds,
          max_depth = models$metadata$hyperparameters$max_depth,
          learning_rate = models$metadata$hyperparameters$learning_rate,
          min_split_loss = models$metadata$hyperparameters$min_split_loss,
          colsample_bytree = models$metadata$hyperparameters$colsample_bytree,
          min_child_weight = models$metadata$hyperparameters$min_child_weight,
          subsample = models$metadata$hyperparameters$subsample,
          objective = models$metadata$hyperparameters$objective,
          partition = NULL, verbose = FALSE
        )
      },
      "net" = {
        fit_abund_net(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          size = models$metadata$hyperparameters$size,
          decay = models$metadata$hyperparameters$decay,
          partition = NULL, verbose = FALSE
        )
      },
      "dnn" = {
        fit_abund_dnn(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          learning_rate = models$metadata$hyperparameters$learning_rate,
          n_epochs = models$metadata$hyperparameters$n_epochs,
          batch_size = models$metadata$hyperparameters$batch_size,
          optimizer = models$metadata$hyperparameters$optimizer,
          validation_patience = models$metadata$hyperparameters$validation_patience,
          fitting_patience = models$metadata$hyperparameters$fitting_patience,
          loss_function = models$metadata$hyperparameters$loss_function,
          weight_decay = models$metadata$hyperparameters$weight_decay,
          custom_architecture = models$metadata$network,
          partition = NULL, verbose = FALSE
        )
      },
      "cnn" = {
        fit_abund_cnn(
          data = db, response = response, predictors = pr_c, predictors_f = pr_f,
          x = extra_args$x, y = extra_args$y, rasters = extra_args$rasters,
          sample_size = extra_args$sample_size,
          learning_rate = if (!is.null(models$optimal_combination$learning_rate)) models$optimal_combination$learning_rate else 0.01,
          n_epochs = if (!is.null(models$optimal_combination$n_epochs)) models$optimal_combination$n_epochs else 10,
          batch_size = if (!is.null(models$optimal_combination$batch_size)) models$optimal_combination$batch_size else 32,
          custom_architecture = extra_args$custom_architecture,
          partition = NULL, verbose = FALSE
        )
      }
    )

    # Predict
    models$model <- m_refit$model
    suppressMessages(p <- adm_predict(models = models, pred = pred, training_data = db))

    # Return cell values as vector
    r_list[[ii]] <- as.vector(p[[1]])
  }
  # parallel::stopCluster(my_cluster)

  # Identify successful iterations
  r_list_clean <- Filter(function(x) is.numeric(x), r_list)
  if (length(r_list_clean) == 0) stop("All bootstrap iterations failed.")

  # Calculate cell-wise standard deviation
  res_mat <- do.call(cbind, r_list_clean)
  unc_vals <- apply(res_mat, 1, stats::sd, na.rm = TRUE)

  # Map back to SpatRaster
  res_raster <- pred[[1]]
  res_raster[] <- unc_vals

  names(res_raster) <- "uncertainty"
  return(res_raster)
}

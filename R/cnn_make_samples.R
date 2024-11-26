#' Creates sample data for Convolutional Neural Network
#'
#' @description This function creates an array of input images and associated responses
#' which can be utilized to train a Convolutional Neural Network.
#'
#' @param data data.frame or tibble. Database that includes longitude, latitude, and response columns.
#' @param x string. Specifying the name of the column with longitude data.
#' @param y string. Specifying the name of the column with latitude data.
#' @param response string. Specifying the name of the column with response.
#' @param raster SpatRaster. Raster from which predictor data will be cropped.
#' @param size numeric. Size of the cropped raster, number o cell in each direction of a focal cell
#' @param raster_padding logical. If TRUE, the raster will be padded when cropping extends beyond its boundaries. Useful for ensuring all focal cells have the same size output even at the edges of the raster. Default FALSE
#' @param padding_method string or NULL. Method used for padding the raster if raster_padding is TRUE. Options are "mean", "median", "zero". Ignored if raster_padding is FALSE. Default NULL
#'
#' @importFrom terra as.array
#'
#' @return A list with two elements - 'predict' (a list of input images) and 'response' (an of response values).
#' Each element in the 'predictors' list is an array representing a cropped image from the input raster.
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
#' cnn_samples <- cnn_make_samples(
#'   data = some_sp,
#'   x = "x", # x coordinates for each point
#'   y = "y", # y coordinates for each point
#'   response = "ind_ha",
#'   raster = envar[[c("bio12","sand","elevation")]],
#'   size = 5 # how many pixels from point to border?
#' )
#' 
#' length(cnn_samples$predictors) #  938 matrix sets
#' dim(cnn_samples$predictors[[1]]) # three 11x11 channels
#' rast(cnn_samples$predictors[[1]]) %>% plot()
#' 
#' cnn_samples$response[[1]] # linked to a label
#' }
cnn_make_samples <- function(data,
                             x,
                             y,
                             response,
                             raster,
                             raster_padding = FALSE,
                             padding_method = NULL,
                             size = 5) {
  predictors <- list()
  responses <- list()
  for (i in 1:nrow(data)) {
    pred_x <- croppin_hood(
      occ = data[i, ],
      x = x,
      y = y,
      raster = raster,
      raster_padding = raster_padding,
      padding_method = padding_method,
      size = size
    )
    
    pred_x <- terra::as.array(pred_x)

    for (j in 1:dim(pred_x)[3]) {
      ary <- pred_x[, , j]
      ary <- ifelse(is.na(ary), mean(ary, na.rm = TRUE), ary)
      pred_x[, , j] <- ary
    }

    if (!is.null(response)){
      pred_y <- data[[i, response]]
    } else if (is.null(response)){
      pred_y <- "null"
    }
    
    predictors <- append(predictors, list(pred_x))
    responses <- append(responses, pred_y)
  }

  if (!is.null(response)){
    data_list <- list(
      "predictors" = predictors,
      "response" = responses
    )
  } else if (is.null(response)){
    data_list <- list(
      "predictors" = predictors
    )
  }
  
  return(data_list)
}

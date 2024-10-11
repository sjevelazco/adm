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
#'
#' @importFrom terra as.array
#'
#' @return A list with two elements - 'predict' (a list of input images) and 'response' (an of response values).
#' Each element in the 'predictors' list is an array representing a cropped image from the input raster.
#'
#' @export
#'
#' @examples
cnn_make_samples <- function(data,
                             x,
                             y,
                             response,
                             raster,
                             size = 5) {
  data <- data

  predictors <- list()
  responses <- list()
  for (i in 1:nrow(data)) {
    pred_x <- croppin_hood(
      occ = data[i, ],
      x = x,
      y = y,
      raster = raster,
      size = size
    )
    pred_x <- terra::as.array(pred_x)

    for (j in 1:dim(pred_x)[3]) {
      ary <- pred_x[, , j]
      ary <- ifelse(is.na(ary), mean(ary, na.rm = TRUE), ary)
      pred_x[, , j] <- ary
    }

    pred_y <- data[[i, response]]

    predictors <- append(predictors, list(pred_x))
    responses <- append(responses, pred_y)
  }

  data_list <- list(
    "predictors" = predictors,
    "response" = responses
  )

  return(data_list)
}

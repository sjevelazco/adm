#' TODO cnn_make_samples
#'
#' @param df
#' @param x character. Column name with longitude data
#' @param y character. Column name with latitude data
#' @param response
#' @param raster
#' @param size
#'
#' @importFrom terra as.array
#' 
#' @return
#' @export
#'
#' @examples
cnn_make_samples <- function(df,
                             x,
                             y,
                             response,
                             raster,
                             size = 5) {
  data <- df

  predictors <- list()
  responses <- list()
  for (i in 1:nrow(data)) {
    pred_x  <- croppin_hood(
      occ = data[i, ],
      x = x,
      y = y,
      raster = raster,
      size = size
    ) 
    pred_x  <- terra::as.array(pred_x)

    for (j in 1:dim(pred_x )[3]) {
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

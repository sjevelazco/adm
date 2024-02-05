#' cnn_make_samples
#'
#' @param df
#' @param longitude
#' @param latitude
#' @param response
#' @param raster
#' @param size
#'
#' @return
#' @export
#'
#' @examples
cnn_make_samples <- function(df,
                             longitude,
                             latitude,
                             response,
                             raster,
                             size = 5) {
  data <- df

  predictors <- list()
  responses <- list()
  for (i in 1:nrow(data)) {
    x <- croppin_hood(
      occ = data[i, ],
      longitude = longitude,
      latitude = latitude,
      raster = raster,
      size = size
    ) %>%
      as.array()

    for (j in 1:dim(x)[3]) {
      ary <- x[, , j]
      ary <- ifelse(is.na(ary), mean(ary, na.rm = TRUE), ary)
      x[, , j] <- ary
    }

    y <- data[[i, response]]

    predictors <- append(predictors, list(x))
    responses <- append(responses, y)
  }

  data_list <- list(
    "predictors" = predictors,
    "response" = responses
  )

  return(data_list)
}

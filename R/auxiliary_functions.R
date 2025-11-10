## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##

#' adapt_df
#'
#' @noRd
adapt_df <- function(data, predictors, predictors_f, response, partition, xy = NULL) {
  data <- data.frame(data)
  if (is.vector(xy)) {
    xy_cols <- data %>%
      dplyr::select(dplyr::all_of(xy))
    xy_cols <- data.frame(xy_cols)
  }

  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::starts_with(partition))
    data <- data.frame(data)
  } else {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), dplyr::starts_with(partition))
    data <- data.frame(data)
    for (i in predictors_f) {
      data[, i] <- as.factor(data[, i])
    }
  }

  if (is.vector(xy)) {
    data <- dplyr::bind_cols(data, xy_cols)
  }

  return(data)
}



#' cnn_get_crop_size
#'
#' @noRd
cnn_get_crop_size <- function(sample_size) {
  if (!is.vector(sample_size)) {
    stop("Please, provide a vector containing the sample size c(width,height)")
  } else if (!(sample_size[[1]] == sample_size[[2]])) {
    stop("adm currently only accepts square samples.")
  } else {
    crop_size <- floor(sample_size[[1]] / 2)
  }

  return(crop_size)
}

#' Crop rasters around a point (Convolutional Neural Networks)
#'
#' @description Crop rasters for a single spatial point. Function used internally to construct Convolutional Neural Networks
#'
#' @param occ tibble or data.frame. Database with response, predictors, and partition values
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param raster SpatRaster. Raster with environmental variables.
#' @param raster_padding logical. If TRUE, the raster will be padded when cropping extends beyond its boundaries. Useful for ensuring all focal cells have the same size output even at the edges of the raster. Default FALSE
#' @param padding_method string or NULL. Method used for padding the raster if raster_padding is TRUE. Options are "mean", "median", "zero". Ignored if raster_padding is FALSE. Default NULL
#' @param size numeric. Size of the cropped raster, number o cell in each direction of a focal cell
#'
#' @importFrom terra colFromX rowFromY xFromCol yFromRow rast ext crop extend
#'
#' @return SpatRaster. Croped raster
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#'
#' # Datasbase with species abundance and x and y coordinates
#' data("sppabund")
#'
#' # Extract data for a single species
#' some_sp <- sppabund %>%
#'   filter(species == "Species three")
#'
#' # Raster data with environmental variables
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)
#'
#' #
#' sampl_r <- croppin_hood(occ = some_sp[1, ], x = "x", y = "y", raster = envar, size = 5)
#' plot(sampl_r)
#' plot(sampl_r[[1]])
#' points(some_sp[1, c("x", "y")], pch = 19)
#' }
croppin_hood <- function(occ,
                         x,
                         y,
                         raster,
                         size,
                         raster_padding = FALSE,
                         padding_method = NULL) {
  if (raster_padding & is.null(padding_method)) {
    stop("Padding method needed.")
  }

  long <- as.numeric(occ[, x])
  lat <- as.numeric(occ[, y])

  if (raster_padding) {
    raster <- terra::extend(raster, c(size, size))
  }

  rst.col <- terra::colFromX(raster, long)
  rst.row <- terra::rowFromY(raster, lat)

  x.max <- terra::xFromCol(raster, rst.col + size)
  x.min <- terra::xFromCol(raster, rst.col - size)
  y.max <- terra::yFromRow(raster, rst.row - size)
  y.min <- terra::yFromRow(raster, rst.row + size)

  r <- terra::rast()
  terra::ext(r) <- c(x.min, x.max, y.min, y.max)

  cropped <- terra::crop(raster, r, snap = "out")

  if (raster_padding) {
    if (padding_method == "mean") {
      padding_values <- terra::global(cropped, mean, na.rm = TRUE)
    } else if (padding_method == "median") {
      padding_values <- terra::global(cropped, median, na.rm = TRUE)
    } else if (padding_method == "zero") {
      padding_values <- data.frame(
        value = rep(0, length(names(cropped))),
        row.names = names(cropped)
      )
    } else {
      stop("Invalid padding method.")
    }

    for (i in 1:length(names(cropped))) {
      cropped[[i]][is.na(cropped[[i]])] <- padding_values[i, ]
    }
  }

  return(cropped)
}

#' Select probability distributions for GAM and GLM
#'
#' @description
#' Select probability distribution available in gamlss.dist suited for a given response variables (e.g., count, zero-inflated) used to fit GAM and GLM models. See \link[gamlss.dist]{gamlss.family} for more details.
#'
#'
#' @param data data.frame or tibble. Database with species abundance
#' @param response character. Column name with species abundance
#'
#' @importFrom dplyr filter select as_tibble arrange
#' @importFrom utils read.delim
#'
#' @return tibble with family_name, family_call, range, and discrete columns (if family distribution
#' or not discrete)
#' @export
#'
#' @examples
#' \dontrun{
#' data(sppabund)
#'
#' family_selector(data = sppabund, response = "ind_ha")
#' }
#'
family_selector <- function(data, response) {
  . <- discrete <- accepts_negatives <- accepts_zero <-
    accepts_one <- one_restricted <- family_name <- family_call <- NULL
  families_bank <-
    system.file("external/families_bank.txt", package = "adm") %>%
    utils::read.delim(., header = TRUE, quote = "\t")

  if (all(round(data[, response]) == data[, response])) {
    # discrete <- TRUE
    message("Response variable is discrete. Both continuous and discrete families will be tested.")
  } else {
    # discrete <- FALSE
    message("Response variable is continuous and need to be transformed to integer in order to test for discrete families.")
    families_bank <- families_bank %>%
      filter(discrete == 0)
  }

  # response has negative values?
  if (any(data[, response] < 0)) {
    families_bank <- families_bank %>%
      filter(accepts_negatives == 1)
  }

  # response haz zeros?
  if (any(data[, response] == 0)) {
    families_bank <- families_bank %>%
      dplyr::filter(accepts_zero == 1)
  }

  # response has ones?
  if (any(data[, response] == 1)) {
    families_bank <- families_bank %>%
      dplyr::filter(accepts_one == 1)
  }

  # response has value > 1?
  if (any(data[, response] > 1)) {
    families_bank <- families_bank %>%
      dplyr::filter(one_restricted == 0)
  }

  testing_families <- families_bank %>%
    dplyr::select(family_name, family_call, range, discrete)

  message("Selected ", nrow(testing_families), " suitable families for the data.")

  return(dplyr::as_tibble(testing_families) %>%
    dplyr::arrange(family_name))
}

#' Calculate the output resolution of a layer
#'
#' @description Calculate the output resolution of a layer or pooling operation in a
#' Convolutional Neural Network.
#'
#' @param type string. Accepted values are "layer" and "pooling".
#' @param in_res integer. It represents the resolution of the input layer.
#' @param kernel_size integer. It refers to the size of the kernel used in the convolution
#' or pooling operation.
#' @param stride integer. It is the stride length for the convolution or pooling operation.
#' Only used when type is "layer".
#' @param padding integer. It is the amount padding added to the input layer.
#' Only used when type is "layer"
#'
#' @return The function returns integer which is the output resolution.
#'
#' @details
#' \itemize{
#' \item When type is "layer" the output resolution is calculated as
#' ((in_res - kernel_size + (2 * padding)) / stride) + 1.
#' \item When type is "pooling", the output resolution is calculated as the floor division
#' of in_res by kernel_size.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Calculating output resolution for a convolution layer
#' res_calculate(type = "layer", in_res = 12, kernel_size = 2, stride = 2, padding = 0)
#'
#' # Calculating output resolution for a pooling layer
#' res_calculate(type = "pooling", in_res = 12, kernel_size = 2)
#' }
res_calculate <-
  function(type = c("layer", "pooling"),
           in_res,
           kernel_size,
           stride,
           padding) {
    type <- match.arg(type)
    if (type == "layer") {
      out_res <- (((in_res - kernel_size + (2 * padding)) / stride) + 1)
    } else if (type == "pooling") {
      out_res <- floor(in_res / kernel_size)
    }

    return(out_res)
  }

#' Construct CNN samples list to use with tune_abund_cnn and fit_abund_cnn
#'
#' @param data 
#' @param x 
#' @param y 
#' @param response 
#' @param folds 
#' @param partition 
#' @param rasters 
#' @param crop_size 
#'
#' @returns a list of arrays
#' @export
#'
#' @examples #TODO
get_partition_samples <- function(data,x,y,response,folds,partition,rasters,crop_size){
  samples_list <- list()
  for (fold in folds) {
    fold_mtx <- data[data[, partition] == fold, c(x, y, response)] %>%
      cnn_make_samples(x, y, response, rasters, size = crop_size) %>%
      list()
    
    names(fold_mtx) <- fold
    
    samples_list <- append(samples_list, fold_mtx)
  }
  
  return(samples_list)
}

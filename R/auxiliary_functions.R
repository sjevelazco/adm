## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##

#' adapt_df
#'
#' @noRd
adapt_df <- function(data, predictors, predictors_f, response, partition, xy = NULL){
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
    data <- bind_cols(data,xy_cols)
  }
  
  return(data)
}


#' pre_tr_te
#'
#' @noRd
pre_tr_te <- function(data, p_names, h) {
  train <- list()
  test <- list()
  
  if (any(c("train", "train-test", "test")
          %in%
          unique(data[, p_names[h]]))) {
    np2 <- 1
    
    filt <- grepl("train", data[, p_names[h]])
    train[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])
    
    filt <- grepl("test", data[, p_names[h]])
    test[[1]] <- data[filt, ] %>%
      dplyr::select(-p_names[!p_names == p_names[h]])
  } else {
    np2 <- max(data[p_names[h]])
    
    for (i in 1:np2) {
      train[[i]] <- data[data[p_names[h]] != i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
      
      test[[i]] <- data[data[p_names[h]] == i, ] %>%
        dplyr::select(-p_names[!p_names == p_names[h]])
    }
  }
  return(list(train = train, test = test, np2 = np2))
}

#' Crop rasters to build samples for Convolutional Neural Networks
#'
#' @param occ
#' @param x
#' @param y
#' @param raster
#' @param size
#'
#' @importFrom terra colFromX rowFromY xFromCol yFromRow rast ext crop
#'
#' @return
#' @export
#'
#' @examples
croppin_hood <- function(occ, x, y, raster, size) {
  long <- as.numeric(occ[, x])
  lat <- as.numeric(occ[, y])

  rst.col <- terra::colFromX(raster, long)
  rst.row <- terra::rowFromY(raster, lat)

  x.max <- terra::xFromCol(raster, rst.col + size)
  x.min <- terra::xFromCol(raster, rst.col - size)
  y.max <- terra::yFromRow(raster, rst.row - size)
  y.min <- terra::yFromRow(raster, rst.row + size)

  r <- terra::rast()
  terra::ext(r) <- c(x.min, x.max, y.min, y.max)

  cropped <- terra::crop(raster, r, snap = "out")
  return(cropped)
}

#' Select probability distributions for GAM and GLM
#'
#' @description
#' Select family distribution suited for a given response variables (e.g., count, zero-inflated) used to fit GAM and GLM models
#'
#'
#' @param data
#' @param response
#'
#' @importFrom dplyr filter select
#' @importFrom utils read.delim
#'
#' @return
#' @export
#'
#' @examples
family_selector <- function(data, response) {
  families_bank <-
    system.file("external/families_bank.txt", package = "adm") %>%
    utils::read.delim(., header = TRUE, quote = "\t")
  
  if (all(round(data[, response]) == data[, response])) {
    # discrete <- TRUE
    message("Response variable is discrete. Both continuous and discrete families will be tested.")
  } else {
    # discrete <- FALSE
    message("Response variable is continuous and need to be rounded in order to test for discrete families.")
    families_bank <- families_bank %>%
      filter(discrete == 0)
  }
  
  # response has negative values?
  if (any(data[,response] < 0)){
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
  
  return(testing_families)
}

#' Select probability distributions for GAM and GLM
#'
#' @description
#'
#'
#'
#' @param type
#' @param in_res
#' @param kernel_size description
#' @param stride description
#' @param padding description
#'
#' @return
#' @export
#'
#' @examples
res_calculate <-
  function(type = c("layer", "pooling"),
           in_res,
           kernel_size,
           stride,
           padding) {
    if (type == "layer") {
      out_res <- (((in_res - kernel_size + (2 * padding)) / stride) + 1)
    } else if (type == "pooling") {
      out_res <- floor(in_res / kernel_size)
    }

    return(out_res)
  }

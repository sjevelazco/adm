## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##


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
    message("Response variable is continuous and will be rounded in order to test for discrete families.")
  }

  # tests if response data has zeros
  if (any(data[, response] == 0)) {
    families_bank <- families_bank %>%
      dplyr::filter(accepts_zero == 1)
  }

  # tests if response data is inside (0,1) interval
  if (any(data[, response] > 1)) {
    families_bank <- families_bank %>%
      dplyr::filter(one_restricted == 0)
  } else if (any(data[, response] == 1)) {
    families_bank <- families_bank %>%
      dplyr::filter(accepts_one == 1)
  } else {
    families_bank <- families_bank %>%
      dplyr::filter(accepts_one == 0)
  }

  testing_families <- families_bank %>%
    dplyr::select(family_call, discrete)

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
           padding)
  {
    if (type == "layer") {
      out_res <- (((in_res - kernel_size + (2 * padding)) / stride) + 1)
    } else if (type == "pooling") {
      out_res <- floor(in_res / kernel_size)
    }
    
    return(out_res)
  }

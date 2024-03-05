## %######################################################%##
#                                                          #
####                Auxiliary functions                 ####
#                                                          #
## %######################################################%##

#' croppin_hood
#' @importFrom terra colFromX rowFromY xFromCol yFromRow rast ext crop
#' @noRd
croppin_hood <- function(occ, longitude, latitude, raster, size) {
  require("terra")

  long <- as.numeric(occ[, longitude])
  lat <- as.numeric(occ[, latitude])

  rst.col <- terra::colFromX(raster, long)
  rst.row <- terra::rowFromY(raster, lat)

  x.max <- terra::xFromCol(raster, rst.col + size)
  x.min <- terra::xFromCol(raster, rst.col - size)
  y.max <- terra::yFromRow(raster, rst.row - size)
  y.min <- terra::yFromRow(raster, rst.row + size)

  r <- terra::rast()
  terra::ext(r) <- c(x.min, x.max, y.min, y.max)

  cropped <- terra::crop(raster, r, snap = "out")
  # plot(cropped)
  # points(x=longitude,y=latitude)
  return(cropped)
}

#' Title
#'
#' @param data 
#' @param response 
#'
#' @return
#' @export
#'
#' @examples
family_selector <-
  function(data,
           response){
    
    data("gamlss_families_table", envir = environment())
    
    if (all(round(data[,response]) == data[,response])) {
      #discrete <- TRUE
      message("Response variable is discrete. Both continuous and discrete families will be tested.")
    } else {
      #discrete <- FALSE
      message("Response variable is continuous and will be rounded in order to test for discrete families.")
    }
    
    # tests if response data has zeros
    if (any(data[,response] == 0)){
      families_bank <- families_bank %>% 
        dplyr::filter(accepts_zero == 1)
    }
    
    # tests if response data is inside (0,1) interval
    if (any(data[,response]>1)){
      families_bank <- families_bank %>% 
        dplyr::filter(one_restricted == 0)
    } else if (any(data[,response]==1)) {
      families_bank <- families_bank %>% 
        dplyr::filter(accepts_one == 1)
    } else {
      families_bank <- families_bank %>% 
        dplyr::filter(accepts_one == 0)
    }
    
    testing_families <- families_bank %>%
      dplyr::select(family_call,discrete)
    
    message("Selected ", nrow(testing_families)," suitable families for the data.")
    
    return(testing_families)
  }

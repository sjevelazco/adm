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

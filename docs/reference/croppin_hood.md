# Crop rasters around a point (Convolutional Neural Networks)

Crop rasters for a single spatial point. Function used internally to
construct Convolutional Neural Networks

## Usage

``` r
croppin_hood(
  occ,
  x,
  y,
  raster,
  size,
  raster_padding = FALSE,
  padding_method = NULL
)
```

## Arguments

- occ:

  tibble or data.frame. Database with response, predictors, and
  partition values

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- raster:

  SpatRaster. Raster with environmental variables.

- size:

  numeric. Size of the cropped raster, number o cell in each direction
  of a focal cell

- raster_padding:

  logical. If TRUE, the raster will be padded when cropping extends
  beyond its boundaries. Useful for ensuring all focal cells have the
  same size output even at the edges of the raster. Default FALSE

- padding_method:

  string or NULL. Method used for padding the raster if raster_padding
  is TRUE. Options are "mean", "median", "zero". Ignored if
  raster_padding is FALSE. Default NULL

## Value

SpatRaster. Croped raster

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)

# Datasbase with species abundance and x and y coordinates
data("sppabund")

# Extract data for a single species
some_sp <- sppabund %>%
  filter(species == "Species three")

# Raster data with environmental variables
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

#
sampl_r <- croppin_hood(occ = some_sp[1, ], x = "x", y = "y", raster = envar, size = 5)
plot(sampl_r)
plot(sampl_r[[1]])
points(some_sp[1, c("x", "y")], pch = 19)
} # }
```

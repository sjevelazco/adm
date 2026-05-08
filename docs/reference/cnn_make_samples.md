# Creates sample data for Convolutional Neural Network

This function creates an array of input images and associated responses
which can be utilized to train a Convolutional Neural Network.

## Usage

``` r
cnn_make_samples(
  data,
  x,
  y,
  response,
  raster,
  raster_padding = FALSE,
  padding_method = NULL,
  size = 5
)
```

## Arguments

- data:

  data.frame or tibble. Database that includes longitude, latitude, and
  response columns.

- x:

  string. Specifying the name of the column with longitude data.

- y:

  string. Specifying the name of the column with latitude data.

- response:

  string. Specifying the name of the column with response.

- raster:

  SpatRaster. Raster from which predictor data will be cropped.

- raster_padding:

  logical. If TRUE, the raster will be padded when cropping extends
  beyond its boundaries. Useful for ensuring all focal cells have the
  same size output even at the edges of the raster. Default FALSE

- padding_method:

  string or NULL. Method used for padding the raster if raster_padding
  is TRUE. Options are "mean", "median", "zero". Ignored if
  raster_padding is FALSE. Default NULL

- size:

  numeric. Size of the cropped raster, number o cell in each direction
  of a focal cell

## Value

A list with two elements - 'predict' (a list of input images) and
'response' (an of response values). Each element in the 'predictors'
list is an array representing a cropped image from the input raster.

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)

# Load data
envar <- system.file("external/envar.tif", package = "adm") %>%
  rast()
data("sppabund")
some_sp <- sppabund %>%
  filter(species == "Species one")

cnn_samples <- cnn_make_samples(
  data = some_sp,
  x = "x", # x coordinates for each point
  y = "y", # y coordinates for each point
  response = "ind_ha",
  raster = envar[[c("bio12", "sand", "elevation")]],
  size = 5 # how many pixels from point to border?
)

length(cnn_samples$predictors) #  938 matrix sets
dim(cnn_samples$predictors[[1]]) # three 11x11 channels
cnn_samples$predictors[[1]] # representing predictor variables
rast(cnn_samples$predictors[[1]]) %>% plot()

cnn_samples$response[[1]] # linked to a label
} # }
```

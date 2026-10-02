# Construct CNN samples list to use with tune_abund_cnn and fit_abund_cnn

Construct CNN samples list to use with tune_abund_cnn and fit_abund_cnn

## Usage

``` r
get_partition_samples(
  data,
  x,
  y,
  response,
  folds,
  partition,
  rasters,
  crop_size
)
```

## Arguments

- data:

  data.frame or tibble. Database with response, coordinates, and
  partition values

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- response:

  character. Column name with species abundance

- folds:

  vector. Values of the partition column (folds) for which samples will
  be built

- partition:

  character. Column name with the partition groups

- rasters:

  SpatRaster. Raster with environmental variables

- crop_size:

  numeric. Number of cells in each direction of a focal cell (see
  [`cnn_make_samples`](https://sjevelazco.github.io/adm/reference/cnn_make_samples.md))

## Value

a list of arrays

# Calculate abundance distribution model uncertainty

This function calculates the uncertainty of an abundance distribution
model by performing a bootstrap procedure. It refits the model multiple
times on resampled data and then calculates the standard deviation of
the predictions across all iterations.

## Usage

``` r
adm_uncertainty(
  models,
  training_data,
  response,
  pred,
  iteration = 50,
  sample_prop = 0.8,
  n_cores = 1
)
```

## Arguments

- models:

  A model object from \`fit_abund\_\*\` or \`tune_abund\_\*\` functions.

- training_data:

  A data.frame or tibble with abundance data and predictors.

- response:

  character. Column name of the response variable.

- pred:

  A SpatRaster object with the environmental layers for projection.

- iteration:

  numeric. The number of bootstrap iterations. Default 50.

- n_cores:

  numeric. The number of cores to use for parallel processing. Default
  1.

- ...:

  Additional arguments passed to refitting functions or
  [`adm_predict`](https://sjevelazco.github.io/adm/reference/adm_predict.md)
  (e.g., `x`, `y`, `rasters`, `sample_size` for CNN;
  `custom_architecture` for DNN/CNN; `invert_transform`,
  `transform_negative` for spatial prediction).

## Value

A SpatRaster object with a single layer representing the model
uncertainty, calculated as the standard deviation of the bootstrap
predictions.

## See also

[`adm_predict`](https://sjevelazco.github.io/adm/reference/adm_predict.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(dplyr)

# Load data
data("cretusa_data")
cretusa_predictors <- system.file("external/cretusa_predictors.tif", package = "adm")
cretusa_predictors <- terra::rast(cretusa_predictors)

species_data <- adm_extract(
  data = cretusa_data,
  x = "x",
  y = "y",
  env_layer = cretusa_predictors
)

# Fit model
mraf <- fit_abund_raf(
  data = species_data,
  response = "ind_ha",
  predictors = c("PC1", "PC2", "PC3"),
  partition = ".part"
)

# Calculate uncertainty
unc <- adm_uncertainty(
  models = mraf,
  training_data = species_data,
  response = "ind_ha",
  pred = cretusa_predictors[[c("PC1", "PC2", "PC3")]],
  iteration = 10,
  n_cores = 2
)

plot(unc)
} # }
```

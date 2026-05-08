# Calculate data to construct bivariate partial dependence plots

Calculate data to construct bivariate partial dependence for two
predictor set

## Usage

``` r
data_abund_bpdp(
  model,
  predictors,
  resolution = 50,
  training_data = NULL,
  invert_transform = NULL,
  response_name = NULL,
  training_boundaries = NULL,
  projection_data = NULL,
  sample_size = NULL,
  training_raster = NULL,
  x_coord = NULL,
  y_coord = NULL
)
```

## Arguments

- model:

  object returned by any fit_abund or tune_abund family functions

- predictors:

  character. Vector with two predictor name(s) to plot. If NULL all
  predictors will be plotted. Default NULL

- resolution:

  numeric. Number of equally spaced points at which to predict
  continuous predictors. Default 50

- training_data:

  data.frame or tibble. Database with response (0,1) and predictor
  values used to fit a model. Default NULL

- invert_transform:

  logical. Invert transformation of response variable. Useful for those
  cases that the response variable was transformed with one of the
  method in
  [`adm_transform`](https://sjevelazco.github.io/adm/reference/adm_transform.md).
  Default NULL

- response_name:

  character. Name of the response variable. Default NULL

- training_boundaries:

  character. Plot training conditions boundaries based on training data
  (i.e., presences, presences and absences, etc). If training_boundaries
  = "convexh", function will delimit training environmental region based
  on a convex-hull. If training_boundaries = "rectangle", function will
  delimit training environmental region based on four straight lines. If
  used any methods it is necessary provide data in training_data
  argument. If NULL all predictors will be used. Default NULL.

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. Default NULL

- sample_size:

  vector. For CNN only. A vector containing the dimensions, in pixels,
  of raster samples. See cnn_make_samples beforehand. Default c(11,11)

- training_raster:

  a terra SpatRaster object. For CNN only. A raster containing the
  predictor variables used in tune_abund_cnn or fit_abund_cnn.

- x_coord:

  character. For CNN only. The name of the column containing longitude
  information for each observation.

- y_coord:

  character. For CNN only. The name of the column containing latitude
  information for each observation.

## Value

A list with two tibbles "pdpdata" and "resid".

- pdpdata: has data to construct partial dependence bivariate plot, the
  first two column includes values of the selected environmental
  variables, the third column the predicted suitability.

- training_boundaries: has data to plot boundaries of training data.

## See also

[`data_abund_pdp`](https://sjevelazco.github.io/adm/reference/data_abund_pdp.md),
[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md),
[`p_abund_bpdp`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)

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

# Fit some models
mglm <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part",
  distribution = "ZAIG",
  poly = 3,
  inter_order = 0,
  predict_part = TRUE
)

# Prepare data for Bivariate Partial Dependence Plots
bpdp_data <- data_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "sand"),
  resolution = 25,
  training_data = some_sp,
  response_name = "Abundance",
  projection_data = envar,
  training_boundaries = "convexh"
)

bpdp_data
} # }
```

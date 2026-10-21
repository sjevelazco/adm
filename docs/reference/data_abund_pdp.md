# Calculate data to construct partial dependence plots

Calculate data to construct partial dependence plots

## Usage

``` r
data_abund_pdp(
  model,
  predictors,
  resolution = 50,
  resid = FALSE,
  training_data = NULL,
  invert_transform = NULL,
  response_name = NULL,
  projection_data = NULL,
  sample_size = c(11, 11),
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

- resid:

  logical. Calculate residuals based on training data. Default FALSE

- training_data:

  data.frame. Database with response and predictor values used to fit a
  model. Default NULL. Required for GLM, GAM, DNN, NET, RAF, SVM models

- invert_transform:

  vector. A vector containing method and terms to invert transformation
  of response variable. Useful for those cases that the response
  variable was transformed with one of the method in
  [`adm_transform`](https://sjevelazco.github.io/adm/reference/adm_transform.md).
  Usage:

  - For "01": invert_transform = c(method = "01", a = min(x), b =
    max(x))

  - For "zscore": invert_transform = c(method = "zscore", a = mean(x), b
    = sd(x))

  - For "log" and "log1: not needed.

  - Can't invert "round" transformations.

  - Default NULL

- response_name:

  character. Name of the response variable. Default NULL

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. When this argument is used, function will calculate
  partial dependence curves distinguishing conditions used in training
  and projection conditions (i.e., projection data present in projection
  area but not training). Default NULL

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

- pdpdata: has data to construct partial dependence plots, the first
  column includes values of the selected environmental variable, the
  second column with predicted suitability, and the third column with
  range type, with two values Training and Projecting, referring to
  suitability calculated within and outside the range of training
  conditions. Third column is only returned if "projection_data"
  argument is used

- resid: has data to plot residuals. The first column includes values of
  the selected environmental variable and the second column with
  predicted suitability.

## See also

[`data_abund_bpdp`](https://sjevelazco.github.io/adm/reference/data_abund_bpdp.md),
[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md),
[`p_abund_bpdp`](https://sjevelazco.github.io/adm/reference/p_abund_bpdp.md)

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)

# Load  data
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
  partition = ".part1",
  distribution = "ZAIG",
  poly = 3,
  inter_order = 0,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco + I(bio12^2) + I(elevation^2) + I(sand^2) + I(bio12^3) + I(elevation^3) + I(sand^3)
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

# Prepare data for Partial Dependence Plots
pdp_data <- data_abund_pdp(
  model = mglm,
  predictors = "bio12",
  resolution = 25,
  resid = TRUE,
  training_data = some_sp,
  response_name = "Abundance",
  projection_data = envar
)

pdp_data
#> $pdpdata
#> # A tibble: 27 × 3
#>      bio12 Abundance Type      
#>      <dbl>     <dbl> <chr>     
#>  1    9.28     0.467 Projection
#>  2  141.       0.956 Projection
#>  3  273.       1.68  Projection
#>  4  405.       2.57  Projection
#>  5  460.       2.97  Training  
#>  6  536.       3.49  Training  
#>  7  668.       4.24  Training  
#>  8  800.       4.68  Training  
#>  9  932.       4.77  Training  
#> 10 1063.       4.54  Training  
#> # ℹ 17 more rows
#> 
#> $resid
#> # A tibble: 938 × 2
#>    bio12 Abundance
#>    <dbl>     <dbl>
#>  1  921.      38.5
#>  2 1008.      40.1
#>  3 1103.      39.3
#>  4 1056.      40.5
#>  5 1175.      38.3
#>  6 2185.      16.2
#>  7 2185.      16.2
#>  8 2403.      15.8
#>  9 2349.      14.6
#> 10 2349.      14.6
#> # ℹ 928 more rows
#> 
# }
```

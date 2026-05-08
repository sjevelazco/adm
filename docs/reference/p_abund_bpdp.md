# Bivariate partial dependence plots for abundance-based distribution models

Create bivariate partial dependence plots to explore the marginal effect
of predictors on modeled abundance

## Usage

``` r
p_abund_bpdp(
  model,
  predictors = NULL,
  resolution = 50,
  training_data = NULL,
  projection_data = NULL,
  training_boundaries = NULL,
  invert_transform = NULL,
  response_name = NULL,
  color_gradient = c("#000004", "#1B0A40", "#4A0C69", "#781B6C", "#A42C5F", "#CD4345",
    "#EC6824", "#FA990B", "#F7CF3D", "#FCFFA4"),
  color_training_boundaries = "white",
  set_max = NULL,
  set_min = NULL,
  theme = ggplot2::theme_classic(),
  sample_size = NULL,
  training_raster = NULL,
  x_coord = NULL,
  y_coord = NULL
)
```

## Arguments

- model:

  A model object found in the first element of the list returned by any
  function from the fit_abund\_ or tune_abund\_ function families

- predictors:

  character. Vector of predictor name(s) to calculate partial dependence
  plots. If NULL all predictors will be used. Default NULL

- resolution:

  numeric. Number of equally spaced points at which to predict abundance
  values for continuous predictors. Default 50

- training_data:

  data.frame or tibble. Database with response and predictor values used
  to fit a model. Default NULL

- projection_data:

  SpatRaster. Raster layer with environmental variables used for model
  projection. When this argument is used, function will calculate
  partial dependence curves distinguishing conditions used in training
  and projection conditions (i.e., projection data present in projection
  area but not training). Default NULL

- training_boundaries:

  character. Plot training conditions boundaries based on training data.
  If training_boundaries = "convexh", function will delimit training
  environmental region based on a convex-hull. If training_boundaries =
  "rectangle", function will delimit training environmental region based
  on four straight lines. If used any methods it is necessary provide
  data in training_data argument.If NULL all predictors will be used.
  Default NULL.

- invert_transform:

  logical. Invert transformation of response variable. Useful for those
  cases that the response variable was transformed with one of the
  method in
  [`adm_transform`](https://sjevelazco.github.io/adm/reference/adm_transform.md).
  Default NULL

- response_name:

  character. Name of the response variable. Default NULL

- color_gradient:

  character. Vector with gradient colors. Default c("#000004",
  "#1B0A40", "#4A0C69", "#781B6C", "#A42C5F", "#CD4345", "#EC6824",
  "#FA990B", "#F7CF3D", "#FCFFA4")

- color_training_boundaries:

  character. A vector with one color used to color points of residuals,
  Default "white"

- set_max:

  numeric. Set a maximum abundance value to plot

- set_min:

  numeric. Set a minimum abundance value to plot

- theme:

  ggplot2 theme. Default ggplot2::theme_classic()

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

A ggplot object

## See also

[`data_abund_pdp`](https://sjevelazco.github.io/adm/reference/data_abund_pdp.md),
[`data_abund_bpdp`](https://sjevelazco.github.io/adm/reference/data_abund_bpdp.md),
[`p_abund_pdp`](https://sjevelazco.github.io/adm/reference/p_abund_pdp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(dplyr)

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

# Bivariate Dependence Plots:
# In different resolutions
p_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "sand"),
  training_data = some_sp,
  resolution = 50
)

p_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "sand"),
  training_data = some_sp,
  resolution = 25
)

# With projection and training boundaries
p_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "elevation", "sand"),
  training_data = some_sp,
  projection_data = envar,
  training_boundaries = "rectangle"
)

p_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "elevation", "sand"),
  training_data = some_sp,
  projection_data = envar,
  training_boundaries = "convexh"
)

# Customize colors and theme
p_abund_bpdp(
  model = mglm,
  predictors = c("bio12", "sand"),
  training_data = some_sp,
  projection_data = envar,
  training_boundaries = "convexh",
  color_gradient =
    c(
      "#122414", "#183C26", "#185437", "#106D43", "#0F874C",
      "#2D9F54", "#61B463", "#8DC982", "#B3E0A7", "#D7F9D0"
    ),
  color_training_boundaries = "purple",
  theme = ggplot2::theme_dark()
)
} # }
```

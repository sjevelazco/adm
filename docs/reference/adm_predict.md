# Spatial predictions from individual and ensemble models

This function allows the geographical prediction of one or more models
constructed with the fit\_ or tune\_ function set, models fitted with
esm\_ function set (i.e., ensemble of small models approach), or models
constructed with fit_ensemble function. It can return continuous or
continuous and binary predictions for one or more thresholds

## Usage

``` r
adm_predict(
  models,
  pred,
  training_data = NULL,
  nchunk = 1,
  predict_area = NULL,
  invert_transform = NULL,
  transform_negative = FALSE,
  sample_size = NULL,
  pred_quantile = 0.5
)
```

## Arguments

- models:

  list of one or more models fitted with fit\_ or tune\_ functions. In
  case use models fitted with fit_ensemble or esm\_ family function only
  one model could be used. Usage models = mglm or models = list(mglm,
  mraf, mgbm)

- pred:

  SpatRaster. Raster layer with predictor variables. Names of layers
  must exactly match those used in model fitting.

- training_data:

  data.frame or tibble. Data used to fit the models. It is necessary to
  predict GAM and GLM models. Default NULL

- nchunk:

  integer. Number of chunks to split data used to predict models (i.e.,
  SpatRaster used in pred argument). Predicting models in chunks helps
  reduce memory requirements in cases where models are predicted for
  large scales and high resolution. Default = 1

- predict_area:

  SpatVector, SpatialPolygon, or SpatialPolygonDataFrame. Spatial
  polygon used for restring prediction into only a given region. Default
  = NULL

- invert_transform:

  named vector. Invert transformation of response variable. Useful for
  those cases that the response variable was transformed with one of the
  method in
  [`adm_transform`](https://sjevelazco.github.io/adm/reference/adm_transform.md).
  Usage = c(method = "anymethod", a = "transformation term a", b =
  "transformation term b"). Default NULL

- transform_negative:

  logical. If TRUE, all negative values in the prediction will be set to
  zero. default FALSE.

- sample_size:

  numeric. A vector containing the dimensions, in pixels, of raster
  samples. See cnn_make_samples beforehand. Default c(11,11)

- pred_quantile:

  numeric. Quantile to predict when `models` is a quantile random forest
  (qrf). Default 0.5

## Value

A list of SpatRaster with continuous and/or binary predictions

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)

data("sppabund")
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Species abundance data, coordinates, and partition
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y, .part1)

# Extract environmental data
some_sp <-
  adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar
  )

# Fit some models
mgam <- fit_abund_gam(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part1",
  distribution = "NO",
  verbose = FALSE
)
#> Formula used for model fitting:
#> ind_ha ~ pb(elevation) + pb(sand) + pb(bio3) + pb(bio12)
#> Warning: additive.fit convergence not obtained in  30  iterations

mraf <- fit_abund_raf(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part1",
  ntree = 100,
  verbose = FALSE
)

# adm_predict can be used for predict one or more models fitted with fit_ or tune_ functions

# a single model
ind_p <- adm_predict(
  models = mgam,
  pred = envar,
  training_data = some_sp,
  nchunk = 1,
  transform_negative = TRUE
)
#> Predicting an individual model
plot(ind_p$gam)


# a list of models
list_p <- adm_predict(
  models = list(mgam, mraf),
  pred = envar,
  training_data = some_sp,
  nchunk = 1,
  transform_negative = TRUE
)
#> Predicting a list of models
list_p
#> $gam
#> class       : SpatRaster
#> size        : 201, 274, 1  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> varname     : envar
#> name        :       gam
#> min value   :         0
#> max value   : 44.503344
#> 
#> $raf
#> class       : SpatRaster
#> size        : 201, 274, 1  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> source(s)   : memory
#> varname     : envar
#> name        :        raf
#> min value   :          0
#> max value   : 148.743333
#> 

# Predicting models in chunks helps reduce memory requirements in
# cases where models are predicted for large scales and high resolution
ind_p <- adm_predict(
  models = mraf,
  pred = envar,
  training_data = some_sp,
  nchunk = 4,
  transform_negative = TRUE
)
#> Predicting an individual model
plot(ind_p$raf)

# }
```

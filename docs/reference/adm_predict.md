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

## Value

A list of SpatRaster with continuous and/or binary predictions

## Examples

``` r
if (FALSE) { # \dontrun{
require(dplyr)
require(terra)

data("sppabund")
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Extract data
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y)

some_sp

some_sp <-
  adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar
  )

# Partition
some_sp <- flexsdm::part_random(
  data = some_sp,
  pr_ab = "ind_ha",
  method = c(method = "rep_kfold", folds = 3, replicates = 3)
)


## %######################################################%##
#                                                          #
####          Create different type of models           ####
#                                                          #
## %######################################################%##
# Fit some models
# require(gamlss)
# m1 <- gamlss::fitDist(some_sp$ind_ha, type="realline")
# m1$fits
# m1$failed
#
# m1 <- gamlss(ind_ha ~ pb(elevation) + pb(sand) + pb(bio3) + pb(bio12), family=NO, data=some_sp)
# choosen_dist <- gamlss::chooseDist(m1, parallel="multicore", ncpus=4, type="realAll")

mgam <- fit_abund_gam(
  data = some_sp %>% dplyr::mutate(ind_ha = round(ind_ha)),
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  predictors_f = "eco",
  partition = ".part",
  distribution = gamlss.dist::PO()
)

mraf <- fit_abund_raf(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part",
)

mgbm <- fit_abund_gbm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part",
  distribution =
  )


## %######################################################%##
#                                                          #
####                        Predict models               ####
#                                                          #
## %######################################################%##

# adm_predict can be used for predict one or more models fitted with fit_ or tune_ functions

# a single model
ind_p <- adm_predict(
  models = mgam,
  pred = envar,
  training_data = some_sp,
  nchunk = 1,
  predict_area = NULL,
  invert_transform = NULL,
  transform_negative = FALSE,
  sample_size = NULL
)

# a list of models
list_p <- adm_predict(
  models = list(mgam, mraf, mgbm),
  pred = envar,
  training_data = some_sp,
  nchunk = 1,
  predict_area = NULL,
  invert_transform = NULL,
  transform_negative = FALSE,
  sample_size = NULL
)

## %######################################################%##
#                                                          #
####              Predict model using chunks            ####
#                                                          #
## %######################################################%##
# Predicting models in chunks helps reduce memory requirements in
# cases where models are predicted for large scales and high resolution

ind_p <- adm_predict(
  models = mraf,
  pred = envar,
  training_data = some_sp,
  nchunk = 4,
  predict_area = NULL,
  invert_transform = NULL,
  transform_negative = FALSE,
)
ind_p
ind_p$raf
plot(ind_p$raf)
} # }
```

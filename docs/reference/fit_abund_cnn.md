# Fit and validate Convolutional Neural Network Model

This function is used to fit a convolutional neural network (CNN) model
for abundance.

## Usage

``` r
fit_abund_cnn(
  data,
  response,
  predictors,
  predictors_f = NULL,
  x,
  y,
  rasters = NULL,
  sample_size,
  partition,
  hold_out_set = NULL,
  predict_part = FALSE,
  learning_rate = 0.01,
  weight_decay = 0,
  n_epochs = 10,
  batch_size = 32,
  validation_patience = 2,
  fitting_patience = 5,
  optimizer = torch::optim_adamw,
  loss_function = torch::nn_l1_loss,
  custom_architecture = NULL,
  samples_list = NULL,
  verbose = TRUE
)
```

## Arguments

- data:

  tibble or data.frame. Database with response, predictors, and
  partition values

- response:

  character. Column name with species abundance.

- predictors:

  character. Vector with the column names of quantitative predictor
  variables (i.e. continuous variables). Usage predictors = c("temp",
  "precipt", "sand")

- predictors_f:

  character. Vector with the column names of qualitative predictor
  variables (i.e. ordinal or nominal variables type). Usage predictors_f
  = c("landform")

- x:

  character. The name of the column containing longitude information for
  each observation.

- y:

  character. The name of the column containing latitude information for
  each observation.

- rasters:

  a terra SpatRaster object. A raster containing the predictor variables
  to be cropped around each observation.

- sample_size:

  numeric. A vector containing the dimensions, in pixels, of raster
  samples. See cnn_make_samples beforehand. Default c(11,11)

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- learning_rate:

  numeric. The size of the step taken during the optimization process.
  Default = 0.01

- n_epochs:

  numeric. Maximum number of times the learning algorithm will work
  through the training set. Default = 10

- batch_size:

  numeric. A batch is a subset of the training set used in a single
  iteration of the training process. The size of each batch is referred
  to as the batch size. Default = 32

- validation_patience:

  numerical. An integer indicating the number of epochs without loss
  improvement tolerated by the algorithm in the validation process. If
  the patience limit is exceeded, the training ends. Default 2

- fitting_patience:

  numerical. The same as validation_patience, but in the final model
  fitting process. Default 5

- custom_architecture:

  a Torch nn_module_generator object. A neural network architecture to
  be used instead of the internal default one. Default NULL

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "luz_module_fitted" object from luz (torch framework). This
  object can be used to predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: Averaged performance metrics (see
  [`adm_eval`](https://sjevelazco.github.io/adm/reference/adm_eval.md)).

- performance_part: Performance metrics for each replica and partition.

- predicted_part: Observed and predicted abundance for each test
  partition.

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Extract data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore reponse variables
some_sp$ind_ha %>% range()
some_sp$ind_ha %>% hist()

# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Generate an architecture
cnn_arch <- generate_cnn_architecture(
  number_of_features = 3,
  number_of_outputs = 1,
  sample_size = c(11, 11),
  number_of_conv_layers = 2,
  conv_layers_size = c(14, 28),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  number_of_fc_layers = 1,
  fc_layers_size = c(28),
  pooling = NULL,
  batch_norm = TRUE,
  dropout = 0,
  verbose = T
)

# Fit a CNN model
mcnn <- fit_abund_cnn(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = NULL,
  partition = ".part",
  x = "x",
  y = "y",
  rasters = envar,
  sample_size = c(11, 11),
  learning_rate = 0.01,
  n_epochs = 100,
  batch_size = 32,
  validation_patience = 2,
  fitting_patience = 5,
  custom_architecture = cnn_arch,
  verbose = TRUE,
  predict_part = TRUE
)

mcnn
} # }
```

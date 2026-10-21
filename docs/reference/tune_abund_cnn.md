# Fit and validate Convolutional Neural Network with exploration of hyper-parameters that optimize performance

Fit and validate Convolutional Neural Network with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_cnn(
  data,
  response,
  predictors,
  predictors_f = NULL,
  x,
  y,
  rasters = NULL,
  samples_list = NULL,
  sample_size,
  partition,
  predict_part = FALSE,
  grid = NULL,
  architectures = NULL,
  metrics = NULL,
  n_cores = 1,
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

  character. Column name with longitude data.

- y:

  character. Column name with latitude data.

- rasters:

  character. Path to the raster file of environmental variables.

- samples_list:

  list. Optional list of pre-computed samples (see
  [`get_partition_samples`](https://sjevelazco.github.io/adm/reference/get_partition_samples.md)),
  used instead of extracting them from `rasters`. Default NULL

- sample_size:

  numeric. The dimension, in pixels, of raster samples. See
  cnn_make_samples beforehand. Default c(11,11)

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- grid:

  tibble or data.frame. A dataframe with "batch_size", "n_epochs",
  "learning_rate" as columns and its values combinations as rows. If no
  grid is provided, function will create a default grid combining the
  next hyperparameters: batch_size = 2^seq(4, 6) n_epochs = 10
  learning_rate = seq(from = 0.1, to = 0.2, by = 0.1)
  validation_patience = 2 fitting_patience = 5. In case one or more
  hyperparameters are provided, the function will complete the grid with
  the default values.

- architectures:

  list or character. A list object containing a list of architectures
  (nn_modules_generators from torch), called "arch_list", and a list of
  matrices describing each architecture, called ("arch_dict"); use
  generate_arch_list function to create it. It's also possible to use
  "fit_intern", what will construct the default neural network
  architecture of fit_abund_cnn. If NULL, a list of architectures will
  be generated. Default NULL

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "luz_module_fitted" object from luz (torch framework). This
  object can be used to predicting.

- predictors: A tibble with quantitative (c column names) and
  qualitative (f column names) variables use for modeling.

- performance: A tibble with selected model's performance metrics
  calculated in adm_eval.

- performance_part: A tibble with performance metrics for each test
  partition.

- predicted_part: A tibble with predicted abundance for each test
  partition.

- optimal_combination: A tibble with the selected hyperparameter
  combination and its performance.

- all_combinations: A tibble with all hyperparameters combinations and
  its performance.

- selected_arch: A numeric vector describing the selected architecture
  layers.

## Examples

``` r
# \donttest{
if (torch::torch_is_installed()) {
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Select data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore response variables
some_sp$ind_ha %>% range()
some_sp$ind_ha %>% hist()

# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# A small subset keeps this example fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 120)), ]

# Generate some architectures
many_archs <- generate_arch_list(
  type = "cnn",
  number_of_features = 3,
  number_of_outputs = 1,
  n_layers = 2,
  n_neurons = c(4, 8),
  sample_size = c(5, 5),
  number_of_fc_layers = c(1), # fully connected layers
  fc_layers_size = c(8),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  batch_norm = FALSE
) %>% select_arch_list(
  type = c("cnn"),
  method = "percentile",
  n_samples = 1,
  min_max = TRUE
)

# Create a grid
# Obs.: the grid is tested with every architecture, thus it can get very large.
cnn_grid <- expand.grid(
  learning_rate = 0.001,
  n_epochs = 3,
  batch_size = c(32),
  validation_patience = 2,
  fitting_patience = 2,
  stringsAsFactors = FALSE
)

# Tune a cnn model
tuned_cnn <- tune_abund_cnn(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = cnn_grid,
  rasters = system.file("external/envar.tif", package = "adm"),
  x = "x",
  y = "y",
  sample_size = c(5, 5),
  architectures = many_archs,
  n_cores = 1
)
tuned_cnn

# It is also possible to use a only one architecture
one_arch <- generate_cnn_architecture(
  number_of_features = 3,
  number_of_outputs = 1,
  sample_size = c(5, 5),
  number_of_conv_layers = 2,
  conv_layers_size = c(4, 8),
  conv_layers_kernel = 3,
  conv_layers_stride = 1,
  conv_layers_padding = 0,
  number_of_fc_layers = 1,
  fc_layers_size = c(8),
  pooling = NULL,
  batch_norm = FALSE,
  dropout = 0,
  verbose = FALSE
)

}

#> Using provided architectures.
#> Adding default hyperparameter for: weight_decay
#> Testing 12 combinations.
#> Searching for optimal hyperparameters...
#> 
#> Fitting the best model...
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> Warning: Some torch operators might not yet be implemented for the MPS device. A
#> temporary fix is to set the `PYTORCH_ENABLE_MPS_FALLBACK=1` to use the CPU as a
#> fall back for those operators:
#> ℹ Add `PYTORCH_ENABLE_MPS_FALLBACK=1` to your `.Renviron` file, for example use
#>   `usethis::edit_r_environ()`.
#> ✖ Using `Sys.setenv()` doesn't work because the env var must be set before R
#>   starts.
#> The best model was achieved with: 
#>  learning_rate = 0.001
#>  patience = 2 and 2, n_epochs = 3, batch_size = 32
# }
```

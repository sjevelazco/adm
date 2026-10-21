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

- hold_out_set:

  tibble or data.frame. Optional hold-out dataset, never used during
  training, with the response and predictor columns. If provided,
  performance metrics are also calculated on it. Default NULL

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- learning_rate:

  numeric. The size of the step taken during the optimization process.
  Default = 0.01

- weight_decay:

  numeric. Weight decay (L2 penalty) used by the optimizer. Default 0

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

- optimizer:

  a torch optimizer generator. Default torch::optim_adamw

- loss_function:

  a torch loss function generator. Default torch::nn_l1_loss

- custom_architecture:

  a Torch nn_module_generator object. A neural network architecture to
  be used instead of the internal default one. Default NULL

- samples_list:

  list. Optional list of pre-computed samples (see
  [`get_partition_samples`](https://sjevelazco.github.io/adm/reference/get_partition_samples.md)),
  used instead of extracting them from `rasters`. Default NULL

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
# \donttest{
if (torch::torch_is_installed()) {
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

# A small subset keeps this example fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 150)), ]

envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Generate an architecture
cnn_arch <- generate_cnn_architecture(
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

# Fit a CNN model
mcnn <- fit_abund_cnn(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = NULL,
  partition = ".part1",
  x = "x",
  y = "y",
  rasters = envar,
  sample_size = c(5, 5),
  learning_rate = 0.01,
  n_epochs = 3,
  batch_size = 32,
  validation_patience = 2,
  fitting_patience = 5,
  custom_architecture = cnn_arch,
  verbose = FALSE,
  predict_part = TRUE
)

mcnn
}

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
#> $model
#> A `luz_module_fitted`
#> ── Time ────────────────────────────────────────────────────────────────────────
#> • Total time: 260ms
#> • Avg time per training epoch: 86ms
#> 
#> ── Results ─────────────────────────────────────────────────────────────────────
#> Metrics observed in the last epoch.
#> 
#> ℹ Training:
#> loss: 23.4182
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> An `nn_module` containing 489 parameters.
#> 
#> ── Modules ─────────────────────────────────────────────────────────────────────
#> • conv1: <nn_conv2d> #112 parameters
#> • conv2: <nn_conv2d> #296 parameters
#> • linear1: <nn_linear> #72 parameters
#> • output: <nn_linear> #9 parameters
#> 
#> $predictors
#> # A tibble: 1 × 5
#>   model response c1    c2        c3   
#>   <chr> <chr>    <chr> <chr>     <chr>
#> 1 cnn   ind_ha   bio12 elevation sand 
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 cnn       25.1   2.93          0.0454         0.122         0.0122
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter    slope  pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl>    <dbl>  <dbl>
#> 1 1       1         cnn    24.8     0.160  -0.000602 28.9  -0.00630 0.0956
#> 2 1       2         cnn    22.3     0.0596  0.180     8.13  2.68    0.0670
#> 3 1       3         cnn    28.1    -0.0830 -0.142    38.0  -0.639   0.223 
#> 
#> $predicted_part
#> # A tibble: 150 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10        2.85
#>  2 1       1             10        3.23
#>  3 1       1             15        5.68
#>  4 1       1             13.3      9.08
#>  5 1       1             40       10.1 
#>  6 1       1             80        7.91
#>  7 1       1             30        7.24
#>  8 1       1             30       10.7 
#>  9 1       1             40        9.23
#> 10 1       1             90       11.4 
#> # ℹ 140 more rows
#> 
#> $metadata
#> $metadata$plataform
#> [1] "aarch64-apple-darwin23"
#> 
#> $metadata$R_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $metadata$run_date
#> [1] "2026-10-21 13:39:05 PDT"
#> 
#> $metadata$lr
#> [1] 0.01
#> 
#> $metadata$weight_decay
#> [1] 0
#> 
#> $metadata$loss
#> An `nn_module` containing 0 parameters.
#> 
#> $metadata$optimizer
#> <optim_adamw> object generator
#>   Inherits from: <inherit>
#>   Public:
#>     initialize: function (params, lr = 0.001, betas = c(0.9, 0.999), eps = 1e-08, 
#>     loop_fun: function (group, param, g, p) 
#>     step: function (closure = NULL) 
#>     clone: function (deep = FALSE) 
#>   Parent env: <environment: 0x6efb33e00>
#>   Locked objects: FALSE
#>   Locked class: FALSE
#>   Portable: TRUE
#> 
#> $metadata$source_function
#> [1] "fit_abund_cnn"
#> 
#> 
# }
```

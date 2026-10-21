# Fit and validate Deep Neural Network model

Fit and validate Deep Neural Network model

## Usage

``` r
fit_abund_dnn(
  data,
  response,
  predictors,
  predictors_f = NULL,
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
  verbose = TRUE,
  learning_monitor = FALSE
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

  numeric. The regularization strength: 0 means no penalty, while higher
  values (e.g. 0.01) apply stronger shrinkage to the weights during
  training. Default is 0

- n_epochs:

  numeric. Max number of times the learning algorithm will work through
  the training set. Default = 10

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

  a torch_optimizer_generator. The optimizer to be used in model
  fitting. Default is torch::optim_adamw.

- loss_function:

  a torch nn_loss. The loss function to be used in model fitting.
  Default is torch::nn_l1_loss.

- custom_architecture:

  a Torch nn_module_generator object or a generate_dnn_architecture
  output. A neural network architecture to be used instead of the
  internal default one. Default NULL

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

- learning_monitor:

  logical. If TRUE, the function will return a tibble containing loss
  values over training epochs. It is useful to create
  learning/convergence plots.

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

# Generate a architecture
dnn_arch <- generate_dnn_architecture(
  number_of_features = 3,
  number_of_outputs = 1,
  number_of_hidden_layers = 3,
  hidden_layers_size = c(4, 8, 4),
  batch_norm = FALSE
)

# Fit a NET model
mdnn <- fit_abund_dnn(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = NULL,
  partition = ".part1",
  learning_rate = 0.01,
  n_epochs = 3,
  batch_size = 32,
  validation_patience = 2,
  fitting_patience = 5,
  custom_architecture = dnn_arch,
  verbose = FALSE,
  predict_part = TRUE
)

mdnn
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
#> • Total time: 266ms
#> • Avg time per training epoch: 88ms
#> 
#> ── Results ─────────────────────────────────────────────────────────────────────
#> Metrics observed in the last epoch.
#> 
#> ℹ Training:
#> loss: 27.4131
#> 
#> ── Model ───────────────────────────────────────────────────────────────────────
#> An `nn_module` containing 97 parameters.
#> 
#> ── Modules ─────────────────────────────────────────────────────────────────────
#> • linear1: <nn_linear> #16 parameters
#> • linear2: <nn_linear> #40 parameters
#> • linear3: <nn_linear> #36 parameters
#> • output: <nn_linear> #5 parameters
#> 
#> $predictors
#> # A tibble: 1 × 5
#>   model response c1    c2        c3   
#>   <chr> <chr>    <chr> <chr>     <chr>
#> 1 dnn   ind_ha   bio12 elevation sand 
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 dnn       24.9   4.75           0.112         0.181        -0.0560
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter  slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl>  <dbl> <dbl>
#> 1 1       1         dnn    23.6     0.194     0.0516  23.5  0.305 0.169
#> 2 1       2         dnn    21.0     0.237    -0.102   30.5 -0.449 0.228
#> 3 1       3         dnn    30.2    -0.0950   -0.117   54.3 -1.17  0.100
#> 
#> $predicted_part
#> # A tibble: 150 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10       16.6 
#>  2 1       1             10       18.6 
#>  3 1       1             15       16.5 
#>  4 1       1             13.3      3.71
#>  5 1       1             40        8.36
#>  6 1       1             80       10.8 
#>  7 1       1             30       10.8 
#>  8 1       1             30       10.8 
#>  9 1       1             40       10.4 
#> 10 1       1             90       10.5 
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
#> [1] "2026-10-21 13:39:06 PDT"
#> 
#> $metadata$network
#> <neural_net> object generator
#>   Inherits from: <inherit>
#>   Public:
#>     .classes: neural_net nn_module
#>     initialize: function () 
#>     forward: function (x) 
#>     clone: function (deep = FALSE, ..., replace_values = TRUE) 
#>   Private:
#>     .__clone_r6__: function (deep = FALSE) 
#>   Parent env: <environment: 0x6f3af50e0>
#>   Locked objects: FALSE
#>   Locked class: FALSE
#>   Portable: TRUE
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$learning_rate
#> [1] 0.01
#> 
#> $metadata$hyperparameters$weight_decay
#> [1] 0
#> 
#> $metadata$hyperparameters$n_epochs
#> [1] 3
#> 
#> $metadata$hyperparameters$batch_size
#> [1] 32
#> 
#> $metadata$hyperparameters$validation_patience
#> [1] 2
#> 
#> $metadata$hyperparameters$fitting_patience
#> [1] 5
#> 
#> $metadata$hyperparameters$optimizer
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
#> $metadata$hyperparameters$loss_function
#> <nn_l1_loss> object generator
#>   Inherits from: <inherit>
#>   Public:
#>     .classes: nn_l1_loss nn_loss nn_module
#>     initialize: function (reduction = "mean") 
#>     forward: function (input, target) 
#>     clone: function (deep = FALSE, ..., replace_values = TRUE) 
#>   Private:
#>     .__clone_r6__: function (deep = FALSE) 
#>   Parent env: <environment: 0x6ef4fee40>
#>   Locked objects: FALSE
#>   Locked class: FALSE
#>   Portable: TRUE
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_dnn"
#> 
#> 
# }
```

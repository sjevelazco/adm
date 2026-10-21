# Fit and validate Artificial Neural Network models

Fit and validate Artificial Neural Network models

## Usage

``` r
fit_abund_net(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  hold_out_set = NULL,
  predict_part = FALSE,
  size,
  decay = 0,
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

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(abund ~ temp + precipt + sand + landform)). Note that the
  variables used here must be consistent with those used in response,
  predictors, and predictors_f arguments. Default NULL

- partition:

  character. Column name with training and validation partition groups.

- hold_out_set:

  tibble or data.frame. Optional hold-out dataset, never used during
  training, with the response and predictor columns. If provided,
  performance metrics are also calculated on it. Default NULL

- predict_part:

  logical. Save predicted abundance for testing data. Default is FALSE.

- size:

  numerical. The size of the hidden layer.

- decay:

  numerial. Value for weight decay. Default 0.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "ksvm" class object from kernlab package. This object can be
  used for predicting.

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
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Extract data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore reponse variables
some_sp$ind_ha %>% range()
#> [1]   0.0000 285.2941
some_sp$ind_ha %>% hist()


# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# Fit a NET model
mnet <- fit_abund_net(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  size = 32,
  decay = 0.1,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

mnet
#> $model
#> a 28-32-1 network with 961 weights
#> inputs: bio12 elevation sand ecoAD ecoAMF ecoAPAF ecoBMDF ecoBY ecoC ecoCADP ecoCAP ecoCDF ecoCM ecoDC ecoE ecoHC ecoHM ecoHP ecoLM ecoP ecoPFS ecoPS ecoSABM ecoSAS ecoSAY ecoSCMS ecoSDMCF ecoUS 
#> output(s): ind_ha 
#> options were - linear output units  decay=0.1
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 net   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 net       23.8   4.06           0.254         0.197          0.232
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear   inter  slope   pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl>   <dbl>  <dbl>   <dbl>
#> 1 1       1         net    22.1     0.421     0.379     2.66  0.819 0.463  
#> 2 1       2         net    20.7     0.304     0.240    12.2   0.546 0.441  
#> 3 1       3         net    28.4     0.0361    0.0763 -339.   13.6   0.00560
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3      17.0
#>  2 1       1             20        17.0
#>  3 1       1             10        17.0
#>  4 1       1             23.3      17.0
#>  5 1       1             10        17.0
#>  6 1       1             10        17.0
#>  7 1       1             20        17.0
#>  8 1       1             10        17.0
#>  9 1       1             15        17.0
#> 10 1       1             10        17.0
#> # ℹ 552 more rows
#> 
#> $metadata
#> $metadata$plataform
#> [1] "aarch64-apple-darwin23"
#> 
#> $metadata$R_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $metadata$run_date
#> [1] "2026-10-21 13:39:16 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6dc42e820>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$size
#> [1] 32
#> 
#> $metadata$hyperparameters$decay
#> [1] 0.1
#> 
#> $metadata$hyperparameters$rang
#> [1] 0.1
#> 
#> $metadata$hyperparameters$maxit
#> [1] 1000
#> 
#> $metadata$hyperparameters$reltol
#> [1] 1e-05
#> 
#> $metadata$hyperparameters$linout
#> [1] TRUE
#> 
#> $metadata$hyperparameters$trace
#> [1] FALSE
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_net"
#> 
#> 
# }
```

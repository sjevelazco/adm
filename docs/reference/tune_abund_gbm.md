# Fit and validate Generalized Boosted Regression models with exploration of hyper-parameters that optimize performance

Fit and validate Generalized Boosted Regression models with exploration
of hyper-parameters that optimize performance

## Usage

``` r
tune_abund_gbm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  partition,
  predict_part = FALSE,
  grid = NULL,
  distribution,
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

- fit_formula:

  formula. A formula object with response and predictor variables (e.g.
  formula(abund ~ temp + precipt + sand + landform)). Note that the
  variables used here must be consistent with those used in response,
  predictors, and predictors_f arguments. Default NULL

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- grid:

  tibble or data.frame. A dataframe with "n.trees", "interaction.depth",
  "n.minobsinnode" and "shrinkage" as columns and its values as rows. If
  no grid is provided, function will create a default grid combining the
  next hyperparameters: n.trees = c(100, 200, 300), interaction.depth =
  c(1, 2, 3), n.minobsinnode = c(5, 10, 15), shrinkage = seq(0.001, 0.1,
  by = 0.05). In case one or more hyperparameters are provided, the
  function will complete the grid with the default values.

- distribution:

  character. A string specifying the distribution to be used. See
  gbm::gbm documentation for details.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gbm" object from gbm package. This object can be used to
  predicting.

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
require(dplyr)

# Database with species abundance and x and y coordinates
data("sppabund")

# Select data for a single species
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(-.part2, -.part3)

# Explore response variables
some_sp$ind_ha %>% range()
#> [1]   0.0000 285.2941
some_sp$ind_ha %>% hist()


# Here we balance number of absences
some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

# A small subset keeps this example fast
set.seed(1)
some_sp <- some_sp[sort(sample(nrow(some_sp), 200)), ]

# Create a grid
gbm_grid <- expand.grid(
  interaction.depth = c(2, 4),
  n.trees = c(50, 100),
  n.minobsinnode = 5,
  shrinkage = c(0.1, 0.5),
  stringsAsFactors = FALSE
)

tuned_gbm <- tune_abund_gbm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = gbm_grid,
  distribution = "gaussian",
  n_cores = 1
)
#> Using provided grid.
#> Searching for optimal hyperparameters...
#> 
#> Fitting the best model...
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> Using 50 trees...
#> -- Partition number 2/3
#> Using 50 trees...
#> -- Partition number 3/3
#> Using 50 trees...
#> The best model was achieved with: 
#>  n.trees = 50, interaction.depth = 4, n.minobsinnode = 5 and shrinkage = 0.1

tuned_gbm
#> $model
#> gbm::gbm(formula = formula1, distribution = distribution, data = data, 
#>     n.trees = n.trees, interaction.depth = interaction.depth, 
#>     n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, bag.fraction = 0.9)
#> A gradient boosted model with gaussian loss function.
#> 50 iterations were performed.
#> There were 4 predictors of which 4 had non-zero influence.
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 gbm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 gbm       25.2   1.08           0.406         0.103          0.305
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         gbm    24.1      0.426     0.298 11.7  0.413 0.722
#> 2 1       2         gbm    26.3      0.295     0.249 13.3  0.509 0.488
#> 3 1       3         gbm    25.1      0.498     0.369  6.46 1.14  0.324
#> 
#> $predicted_part
#> # A tibble: 200 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10        29.1
#>  2 1       1             10        11.5
#>  3 1       1             10        18.6
#>  4 1       1             15        20.0
#>  5 1       1             13.3      32.1
#>  6 1       1             40        29.5
#>  7 1       1             10        38.5
#>  8 1       1             80        28.1
#>  9 1       1             30        53.4
#> 10 1       1             30        38.8
#> # ℹ 190 more rows
#> 
#> $metadata
#> $metadata$plataform
#> [1] "aarch64-apple-darwin23"
#> 
#> $metadata$R_version
#> [1] "R version 4.6.1 (2026-06-24)"
#> 
#> $metadata$run_date
#> [1] "2026-10-21 13:40:22 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6d01d8708>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "gaussian"
#> 
#> $metadata$hyperparameters$n.trees
#> [1] 50
#> 
#> $metadata$hyperparameters$interaction.depth
#> [1] 4
#> 
#> $metadata$hyperparameters$n.minobsinnode
#> [1] 5
#> 
#> $metadata$hyperparameters$shrinkage
#> [1] 0.1
#> 
#> $metadata$hyperparameters$bag.fraction
#> [1] 0.9
#> 
#> 
#> $metadata$source_function
#> [1] "tune_abund_gbm"
#> 
#> 
#> $optimal_combination
#>   comb_id interaction.depth n.trees n.minobsinnode shrinkage model mae_mean
#> 2  comb_2                 4      50              5       0.1   gbm  24.5396
#>     mae_sd corr_spear_mean corr_spear_sd corr_pear_mean corr_pear_sd inter_mean
#> 2 1.175909        0.439382     0.0942566      0.3188383   0.05439596   8.543514
#>   inter_sd slope_mean  slope_sd pdisp_mean  pdisp_sd
#> 2 4.989447  0.7635458 0.4627265  0.4925684 0.2032781
#> 
#> $all_combinations
#>   comb_id interaction.depth n.trees n.minobsinnode shrinkage model mae_mean
#> 1  comb_1                 2      50              5       0.1   gbm 24.59903
#> 2  comb_2                 4      50              5       0.1   gbm 24.53960
#> 3  comb_3                 2     100              5       0.1   gbm 25.30358
#> 4  comb_4                 4     100              5       0.1   gbm 25.51844
#> 5  comb_5                 2      50              5       0.5   gbm 25.55495
#> 6  comb_6                 4      50              5       0.5   gbm 28.23054
#> 7  comb_7                 2     100              5       0.5   gbm 26.89802
#> 8  comb_8                 4     100              5       0.5   gbm 28.69277
#>      mae_sd corr_spear_mean corr_spear_sd corr_pear_mean corr_pear_sd
#> 1 1.3739845       0.4212312    0.10987698      0.3178379   0.03956656
#> 2 1.1759085       0.4393820    0.09425660      0.3188383   0.05439596
#> 3 1.7393081       0.3899737    0.04867160      0.2945648   0.01674861
#> 4 1.6921530       0.3808268    0.07164847      0.2837482   0.01158143
#> 5 1.6824385       0.3660762    0.06735540      0.2750621   0.01010270
#> 6 2.2534463       0.3006351    0.03545816      0.2539959   0.04217415
#> 7 0.3413461       0.2989221    0.05233153      0.2615484   0.03489487
#> 8 1.3733684       0.3067919    0.03729949      0.2463326   0.02059558
#>   inter_mean  inter_sd slope_mean   slope_sd pdisp_mean  pdisp_sd
#> 1   6.288299 3.3121072  0.8502268 0.43747738  0.4296063 0.1818291
#> 2   8.543514 4.9894471  0.7635458 0.46272649  0.4925684 0.2032781
#> 3  11.108509 0.3570896  0.6610716 0.29073631  0.4995977 0.1950403
#> 4  13.720162 0.4256199  0.5663821 0.25134430  0.5622595 0.2142053
#> 5  16.498663 2.4877867  0.4844660 0.16721275  0.6140074 0.2056581
#> 6  19.680266 3.1836586  0.3580724 0.07271545  0.7350929 0.2302327
#> 7  17.405842 1.4834507  0.4367597 0.19792706  0.6754013 0.2558619
#> 8  19.691298 3.2110553  0.3739251 0.16027813  0.7313334 0.2629399
#> 
# }
```

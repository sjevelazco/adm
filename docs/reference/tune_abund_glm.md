# Fit and validate Generalized Linear Models with exploration of hyper-parameters that optimize performance

Fit and validate Generalized Linear Models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_glm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  sigma_formula = ~1,
  nu_formula = ~1,
  tau_formula = ~1,
  partition,
  predict_part = FALSE,
  grid = NULL,
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

- sigma_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  sigma_formula = ~ precipt + temp

- nu_formula:

  formula. formula for fitting a model to the nu parameter. Usage
  nu_formula = ~ precipt + temp

- tau_formula:

  formula. formula for fitting a model to the tau parameter. Usage
  tau_formula = ~ precipt + temp

- partition:

  character. Column name with training and validation partition groups.

- predict_part:

  logical. Save predicted abundance for testing data. Default = FALSE

- grid:

  tibble or data.frame. A dataframe with "distribution", "poly",
  "inter_order" as columns and its values combinations as rows. If no
  grid is provided, function will create a default grid combining the
  next hyperparameters: poly = c(1, 2, 3), inter_order = c(0, 1, 2),
  distribution = families_hp\$family_call. In case one or more
  hyperparameters are provided, the function will complete the grid with
  the default values.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope").

- n_cores:

  numeric. Number of cores used in parallel processing.

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gamlss" object from gamlss package. This object can be used
  to predicting.

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

## Examples

``` r
# \donttest{
require(dplyr)
require(gamlss)

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
some_sp <- some_sp[sort(sample(nrow(some_sp), 150)), ]

# Explore different family distributions
suitable_distributions <- family_selector(data = some_sp, response = "ind_ha")
#> Response variable is continuous and need to be transformed to integer in order to test for discrete families.
#> Selected 29 suitable families for the data.
suitable_distributions
#> # A tibble: 29 × 4
#>    family_name                         family_call range   discrete
#>    <chr>                               <chr>       <chr>      <int>
#>  1 Exponential Gaussian                exGAUS      (-∞, ∞)        0
#>  2 Exponential generalized Beta type 2 EGB2        (-∞, ∞)        0
#>  3 Generalized t                       GT          (-∞, ∞)        0
#>  4 Gumbel                              GU          (-∞, ∞)        0
#>  5 Johnson's SU                        JSU         (-∞, ∞)        0
#>  6 Logistic                            LO          (-∞, ∞)        0
#>  7 Normal                              NO          (-∞, ∞)        0
#>  8 Normal Exponential t                NET         (-∞, ∞)        0
#>  9 Normal Family                       NOF         (-∞, ∞)        0
#> 10 Normal Linear Quadratic             LQNO        (-∞, ∞)        0
#> # ℹ 19 more rows

# Create a grid
glm_grid <- expand.grid(
  poly = c(1, 2),
  inter_order = 1,
  distribution = c("NO", "TF"),
  stringsAsFactors = FALSE
)

# Tune a GLM model
tuned_glm <- tune_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  fit_formula = formula("ind_ha ~ bio12 + elevation + sand + eco"),
  sigma_formula = formula("ind_ha ~ bio12 + elevation"),
  nu_formula = formula("ind_ha ~ bio12 + elevation"),
  predictors_f = c("eco"),
  partition = ".part1",
  predict_part = TRUE,
  metrics = c("corr_pear", "mae"),
  grid = glm_grid,
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
#> -- Partition number 2/3
#> -- Partition number 3/3
#> The best model was achieved with:
#>  distribution = NO
#>  poly = 1
#>  inter_order = 1

tuned_glm
#> $model
#> 
#> Family:  c("NO", "Normal") 
#> Fitting method: RS() 
#> 
#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  
#>     nu.formula = nu_formula, tau.formula = tau_formula,  
#>     family = family, data = data, control = control_gamlss,      trace = FALSE) 
#> 
#> 
#> Mu Coefficients:
#> (Intercept)        bio12    elevation         sand        ecoAD       ecoAMF  
#>   269.77547     -0.02078     -0.02851     -0.03579           NA   -185.41141  
#>     ecoAPAF      ecoBMDF        ecoBY         ecoC      ecoCADP       ecoCAP  
#>  -195.22335           NA           NA   -192.21362           NA           NA  
#>      ecoCDF        ecoCM        ecoDC         ecoE        ecoHC        ecoHM  
#>          NA           NA   -210.22284           NA   -182.34000   -164.05490  
#>       ecoHP        ecoLM         ecoP       ecoPFS        ecoPS      ecoSABM  
#>          NA           NA           NA           NA           NA           NA  
#>      ecoSAS       ecoSAY      ecoSCMS     ecoSDMCF        ecoUS  
#>          NA   -177.13422           NA           NA           NA  
#> Sigma Coefficients:
#> (Intercept)        bio12    elevation  
#>   4.5414032   -0.0006673   -0.0008272  
#> 
#>  Degrees of Freedom for the fit: 14 Residual Deg. of Freedom   136 
#> Global Deviance:     1420.24 
#>             AIC:     1448.24 
#>             SBC:     1490.39 
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 glm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 glm       25.2   4.50           0.201         0.132          0.255
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         glm    26.4     0.185      0.119 23.2  0.181 0.658
#> 2 1       2         glm    20.2     0.340      0.417  1.06 1.01  0.412
#> 3 1       3         glm    29.0     0.0767     0.230  5.24 0.707 0.325
#> 
#> $predicted_part
#> # A tibble: 150 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10        22.8
#>  2 1       1             10        16.9
#>  3 1       1             15        22.8
#>  4 1       1             13.3      25.5
#>  5 1       1             40        37.3
#>  6 1       1             80        20.7
#>  7 1       1             30        20.1
#>  8 1       1             30        20.1
#>  9 1       1             40        19.6
#> 10 1       1             90        50.1
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
#> [1] "2026-10-21 13:40:26 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6f14ba1c8>
#> 
#> $metadata$sigma.formula
#> ind_ha ~ bio12 + elevation
#> <environment: 0x6f14ba1c8>
#> 
#> $metadata$nu.formula
#> ind_ha ~ bio12 + elevation
#> <environment: 0x6f14ba1c8>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6f47cf770>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "NO"
#> 
#> $metadata$hyperparameters$poly
#> [1] 1
#> 
#> $metadata$hyperparameters$inter_order
#> [1] 1
#> 
#> 
#> $metadata$source_function
#> [1] "tune_abund_glm"
#> 
#> 
#> $optimal_combination
#>   comb_id poly inter_order distribution model mae_mean   mae_sd corr_spear_mean
#> 1  comb_1    1           1           NO   glm 25.21718 4.501824       0.2007036
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1     0.1324206      0.2554389    0.1506896   9.836541 11.77087  0.6341309
#>    slope_sd pdisp_mean  pdisp_sd
#> 1 0.4210576  0.4648498 0.1729726
#> 
#> $all_combinations
#>   comb_id poly inter_order distribution model mae_mean   mae_sd corr_spear_mean
#> 1  comb_1    1           1           NO   glm 25.21718 4.501824       0.2007036
#> 2  comb_2    2           1           NO   glm 25.21718 4.501824       0.2007036
#> 3  comb_3    1           1           TF   glm 24.54533 3.877920       0.1253235
#> 4  comb_4    2           1           TF   glm 24.54533 3.877920       0.1253235
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1     0.1324206      0.2554389    0.1506896   9.836541 11.77087  0.6341309
#> 2     0.1324206      0.2554389    0.1506896   9.836541 11.77087  0.6341309
#> 3     0.1303364      0.1909297    0.2537049  14.973462 14.22892  0.6558709
#> 4     0.1303364      0.1909297    0.2537049  14.973462 14.22892  0.6558709
#>    slope_sd pdisp_mean  pdisp_sd
#> 1 0.4210576  0.4648498 0.1729726
#> 2 0.4210576  0.4648498 0.1729726
#> 3 0.7848827  0.3978921 0.1716193
#> 4 0.7848827  0.3978921 0.1716193
#> 
# }
```

# Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance

Fit and validate Generalized Additive Models with exploration of
hyper-parameters that optimize performance

## Usage

``` r
tune_abund_gam(
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

  tibble or data.frame. A dataframe with 'distribution', 'inter' as
  columns and its values combinations as rows. If no grid is provided,
  function will create a default grid combining the next
  hyperparameters: distribution = all families selected by
  family_selector, inter = "automatic". In case one or more
  hyperparameters are provided, the function will complete the grid with
  the default values.

- metrics:

  character. Vector with one or more metrics from
  c("corr_spear","corr_pear","mae","pdisp","inter","slope")

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
some_sp <- some_sp[sort(sample(nrow(some_sp), 200)), ]
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
gam_grid <- expand.grid(
  distribution = c("NO", "TF"),
  stringsAsFactors = FALSE
)
# Tune a GAM model
tuned_gam <- tune_abund_gam(
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
  grid = gam_grid,
  n_cores = 1
)
#> Adding default hyperparameter for: inter
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
#>  family = NO
#>  inter = automatic

tuned_gam
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
#>   2.380e+02   -2.654e-03   -3.976e-02   -4.495e-02           NA   -1.885e+02  
#>     ecoAPAF      ecoBMDF        ecoBY         ecoC      ecoCADP       ecoCAP  
#>  -1.917e+02           NA           NA   -1.868e+02           NA           NA  
#>      ecoCDF        ecoCM        ecoDC         ecoE        ecoHC        ecoHM  
#>          NA           NA   -1.942e+02           NA   -1.723e+02           NA  
#>       ecoHP        ecoLM         ecoP       ecoPFS        ecoPS      ecoSABM  
#>          NA           NA   -2.148e+02           NA           NA           NA  
#>      ecoSAS       ecoSAY      ecoSCMS     ecoSDMCF        ecoUS  
#>          NA   -1.359e+02           NA           NA           NA  
#> Sigma Coefficients:
#> (Intercept)  
#>       3.495  
#> 
#>  Degrees of Freedom for the fit: 12 Residual Deg. of Freedom   188 
#> Global Deviance:     1965.63 
#>             AIC:     1989.63 
#>             SBC:     2029.21 
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 gam   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1 gam       23.8   1.84           0.430        0.0329          0.359
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         gam    22.1      0.419     0.373  9.13 0.651 0.573
#> 2 1       2         gam    25.8      0.404     0.138 22.5  0.232 0.593
#> 3 1       3         gam    23.5      0.467     0.567 -6.49 1.29  0.438
#> 
#> $predicted_part
#> # A tibble: 200 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             10        24.1
#>  2 1       1             10        11.0
#>  3 1       1             10        18.3
#>  4 1       1             15        25.3
#>  5 1       1             13.3      43.2
#>  6 1       1             40        34.5
#>  7 1       1             10        27.2
#>  8 1       1             80        25.7
#>  9 1       1             30        25.0
#> 10 1       1             30        25.1
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
#> [1] "2026-10-21 13:40:19 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6f00bb620>
#> 
#> $metadata$sigma.formula
#> ~1
#> <environment: 0x6f5384778>
#> 
#> $metadata$nu.formula
#> ~1
#> <environment: 0x6f5384778>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6f5384778>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "NO"
#> 
#> $metadata$hyperparameters$inter
#> [1] "automatic"
#> 
#> 
#> $metadata$source_function
#> [1] "tune_abund_gam"
#> 
#> 
#> $optimal_combination
#>   comb_id     inter distribution model mae_mean   mae_sd corr_spear_mean
#> 1  comb_1 automatic           NO   gam 23.81795 1.838511       0.4300649
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1    0.03290869      0.3591002    0.2149042   8.377114 14.50427  0.7257494
#>    slope_sd pdisp_mean   pdisp_sd
#> 1 0.5354011  0.5347839 0.08465585
#> 
#> $all_combinations
#>   comb_id     inter distribution model mae_mean   mae_sd corr_spear_mean
#> 1  comb_1 automatic           NO   gam 23.81795 1.838511       0.4300649
#> 2  comb_2 automatic           TF   gam 23.16517 1.993356       0.4199715
#>   corr_spear_sd corr_pear_mean corr_pear_sd inter_mean inter_sd slope_mean
#> 1    0.03290869      0.3591002    0.2149042   8.377114 14.50427  0.7257494
#> 2    0.05375432      0.2833611    0.2773289  11.094396 21.29107  1.5590959
#>    slope_sd pdisp_mean   pdisp_sd
#> 1 0.5354011  0.5347839 0.08465585
#> 2 2.0653894  0.3626862 0.20965705
#> 
# }
```

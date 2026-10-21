# Fit and validate Generalized Linear Models

Fit and validate Generalized Linear Models

## Usage

``` r
fit_abund_glm(
  data,
  response,
  predictors,
  predictors_f = NULL,
  fit_formula = NULL,
  sigma_formula = ~1,
  nu_formula = ~1,
  tau_formula = ~1,
  partition = NULL,
  hold_out_set = NULL,
  predict_part = FALSE,
  distribution = NULL,
  poly = 0,
  inter_order = 0,
  control_gamlss = gamlss::gamlss.control(trace = FALSE),
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

- hold_out_set:

  tibble or data.frame. Optional hold-out dataset, never used during
  training, with the response and predictor columns. If provided,
  performance metrics are also calculated on it. Default NULL

- predict_part:

  logical. Save predicted abundance for testing data. Default is FALSE.

- distribution:

  character. A string specifying the distribution to be used. See
  [gamlss.family](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html)
  documentation for details. Use distribution = gamlss.dist::NO().
  Default NULL

- poly:

  integer \>= 2. If used with values \>= 2 model will use polynomials
  for those continuous variables (i.e. used in predictors argument).
  Default is 0.

- inter_order:

  integer \>= 0. The interaction order between explanatory variables.
  Default is 0.

- control_gamlss:

  function. control parameters of the outer iterations algorithm in
  gamlss See
  [gamlss.control](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  documentation for details. Default gamlss.control()

- verbose:

  logical. If FALSE, disables all console messages. Default TRUE

## Value

A list object with:

- model: A "gamlss" class object from gamlss package. This object can be
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
require(terra)
require(dplyr)
require(gamlss)

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

# Explore different family distributions
family_selector(data = some_sp, response = "ind_ha") %>% tail()
#> Response variable is continuous and need to be transformed to integer in order to test for discrete families.
#> Selected 29 suitable families for the data.
#> # A tibble: 6 × 4
#>   family_name         family_call range   discrete
#>   <chr>               <chr>       <chr>      <int>
#> 1 Skew t type 3       ST3         (-∞, ∞)        0
#> 2 Skew t type 4       ST4         (-∞, ∞)        0
#> 3 Skew t type 5       ST5         (-∞, ∞)        0
#> 4 Zero Adjusted Gamma ZAGA        [0,∞)          0
#> 5 Zero adjusted IG    ZAIG        [0,∞)          0
#> 6 t-distribution      TF          (-∞, ∞)        0

# Fit a GLM model
glm_1 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  distribution = "ZAGA",
  poly = 0,
  inter_order = 0,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

glm_1
#> $model
#> 
#> Family:  c("ZAGA", "Zero Adjusted GA") 
#> Fitting method: RS() 
#> 
#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  
#>     nu.formula = nu_formula, tau.formula = tau_formula,  
#>     family = family, data = data, control = control_gamlss,      trace = FALSE) 
#> 
#> 
#> Mu Coefficients:
#> (Intercept)        bio12    elevation         sand        ecoAD       ecoAMF  
#>   5.7846381   -0.0005815   -0.0007766   -0.0004089           NA   -1.3860558  
#>     ecoAPAF      ecoBMDF        ecoBY         ecoC      ecoCADP       ecoCAP  
#>  -1.3524058           NA           NA   -1.2336582           NA    1.4524114  
#>      ecoCDF        ecoCM        ecoDC         ecoE        ecoHC        ecoHM  
#>          NA           NA   -1.4433830   -2.0835995   -0.9277679           NA  
#>       ecoHP        ecoLM         ecoP       ecoPFS        ecoPS      ecoSABM  
#>          NA           NA           NA   -1.8347252           NA           NA  
#>      ecoSAS       ecoSAY      ecoSCMS     ecoSDMCF        ecoUS  
#>          NA   -0.6534990           NA           NA           NA  
#> Sigma Coefficients:
#> (Intercept)  
#>      -0.212  
#> Nu Coefficients:
#> (Intercept)  
#>      -1.605  
#> 
#>  Degrees of Freedom for the fit: 15 Residual Deg. of Freedom   547 
#> Global Deviance:     4586.79 
#>             AIC:     4616.79 
#>             SBC:     4681.76 
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
#> 1 glm       25.6   1.43           0.295         0.116          0.282
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         glm    26.8      0.224     0.202 13.7  0.364 0.554
#> 2 1       2         glm    24.0      0.232     0.242  7.94 0.543 0.446
#> 3 1       3         glm    26.0      0.429     0.403 -8.58 1.28  0.316
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3      36.9
#>  2 1       1             20        14.7
#>  3 1       1             10        11.6
#>  4 1       1             23.3      13.9
#>  5 1       1             10        28.4
#>  6 1       1             10        26.0
#>  7 1       1             20        20.9
#>  8 1       1             10        16.5
#>  9 1       1             15        16.3
#> 10 1       1             10        15.0
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
#> [1] "2026-10-21 13:39:14 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6fb4f8d28>
#> 
#> $metadata$sigma.formula
#> ~1
#> <environment: 0x6fb4f8d28>
#> 
#> $metadata$nu.formula
#> ~1
#> <environment: 0x6fb4f8d28>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6fb4f8d28>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "ZAGA"
#> 
#> $metadata$hyperparameters$poly
#> [1] 0
#> 
#> $metadata$hyperparameters$inter_order
#> [1] 0
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_glm"
#> 
#> 

# Using second order polynomials and first order interaction terms
glm_2 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  distribution = "ZAGA",
  poly = 2,
  inter_order = 1,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco + I(bio12^2) + I(elevation^2) + I(sand^2) + bio12:elevation + bio12:sand + bio12:eco + elevation:sand + eco:elevation + eco:sand
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

glm_2
#> $model
#> 
#> Family:  c("ZAGA", "Zero Adjusted GA") 
#> Fitting method: RS() 
#> 
#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  
#>     nu.formula = nu_formula, tau.formula = tau_formula,  
#>     family = family, data = data, control = control_gamlss,      trace = FALSE) 
#> 
#> 
#> Mu Coefficients:
#>        (Intercept)               bio12           elevation                sand  
#>          1.151e+03          -6.218e-01           8.841e-01          -4.165e-01  
#>              ecoAD              ecoAMF             ecoAPAF             ecoBMDF  
#>                 NA          -1.149e+03          -1.146e+03                  NA  
#>              ecoBY                ecoC             ecoCADP              ecoCAP  
#>                 NA          -1.145e+03                  NA          -2.983e+03  
#>             ecoCDF               ecoCM               ecoDC                ecoE  
#>                 NA                  NA          -1.149e+03          -1.383e+02  
#>              ecoHC               ecoHM               ecoHP               ecoLM  
#>         -1.146e+03                  NA                  NA                  NA  
#>               ecoP              ecoPFS               ecoPS             ecoSABM  
#>                 NA          -1.265e+03                  NA                  NA  
#>             ecoSAS              ecoSAY             ecoSCMS            ecoSDMCF  
#>                 NA          -1.144e+03                  NA                  NA  
#>              ecoUS          I(bio12^2)      I(elevation^2)           I(sand^2)  
#>                 NA          -5.325e-08           6.700e-07           2.455e-06  
#>    bio12:elevation          bio12:sand         bio12:ecoAD        bio12:ecoAMF  
#>          1.263e-06           3.309e-07                  NA           6.221e-01  
#>      bio12:ecoAPAF       bio12:ecoBMDF         bio12:ecoBY          bio12:ecoC  
#>          6.211e-01                  NA                  NA           6.188e-01  
#>      bio12:ecoCADP        bio12:ecoCAP        bio12:ecoCDF         bio12:ecoCM  
#>                 NA                  NA                  NA                  NA  
#>        bio12:ecoDC          bio12:ecoE         bio12:ecoHC         bio12:ecoHM  
#>          6.226e-01                  NA           6.210e-01                  NA  
#>        bio12:ecoHP         bio12:ecoLM          bio12:ecoP        bio12:ecoPFS  
#>                 NA                  NA                  NA           6.993e-01  
#>        bio12:ecoPS       bio12:ecoSABM        bio12:ecoSAS        bio12:ecoSAY  
#>                 NA                  NA                  NA           6.205e-01  
#>      bio12:ecoSCMS      bio12:ecoSDMCF         bio12:ecoUS      elevation:sand  
#>                 NA                  NA                  NA          -8.723e-06  
#>    elevation:ecoAD    elevation:ecoAMF   elevation:ecoAPAF   elevation:ecoBMDF  
#>                 NA          -8.875e-01          -8.893e-01                  NA  
#>    elevation:ecoBY      elevation:ecoC   elevation:ecoCADP    elevation:ecoCAP  
#>                 NA          -8.796e-01                  NA                  NA  
#>   elevation:ecoCDF     elevation:ecoCM     elevation:ecoDC      elevation:ecoE  
#>                 NA                  NA          -8.846e-01                  NA  
#>    elevation:ecoHC     elevation:ecoHM     elevation:ecoHP     elevation:ecoLM  
#>         -8.866e-01                  NA                  NA                  NA  
#>     elevation:ecoP    elevation:ecoPFS     elevation:ecoPS   elevation:ecoSABM  
#>                 NA          -7.889e-01                  NA                  NA  
#>   elevation:ecoSAS    elevation:ecoSAY   elevation:ecoSCMS  elevation:ecoSDMCF  
#>                 NA          -8.842e-01                  NA                  NA  
#>    elevation:ecoUS          sand:ecoAD         sand:ecoAMF        sand:ecoAPAF  
#>                 NA                  NA           4.204e-01           4.162e-01  
#>       sand:ecoBMDF          sand:ecoBY           sand:ecoC        sand:ecoCADP  
#>                 NA                  NA           4.164e-01                  NA  
#>        sand:ecoCAP         sand:ecoCDF          sand:ecoCM          sand:ecoDC  
#>                 NA                  NA                  NA           4.195e-01  
#>          sand:ecoE          sand:ecoHC          sand:ecoHM          sand:ecoHP  
#>                 NA           4.138e-01                  NA                  NA  
#>         sand:ecoLM           sand:ecoP         sand:ecoPFS          sand:ecoPS  
#>                 NA                  NA           4.068e-01                  NA  
#>       sand:ecoSABM         sand:ecoSAS         sand:ecoSAY        sand:ecoSCMS  
#>                 NA                  NA           4.134e-01                  NA  
#>      sand:ecoSDMCF          sand:ecoUS  
#>                 NA                  NA  
#> Sigma Coefficients:
#> (Intercept)  
#>     -0.2606  
#> Nu Coefficients:
#> (Intercept)  
#>      -1.605  
#> 
#>  Degrees of Freedom for the fit: 42 Residual Deg. of Freedom   520 
#> Global Deviance:     4532.49 
#>             AIC:     4616.49 
#>             SBC:     4798.42 
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 glm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean  mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>   <dbl>           <dbl>         <dbl>          <dbl>
#> 1 glm    1.88e18 3.26e18           0.232         0.160         0.0221
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model     mae corr_spear corr_pear inter     slope   pdisp
#>   <chr>   <chr>     <chr>   <dbl>      <dbl>     <dbl> <dbl>     <dbl>   <dbl>
#> 1 1       1         glm   5.65e18     0.0561   -0.0579  27.1 -2.50e-20 2.32e18
#> 2 1       2         glm   2.87e 1     0.274     0.0899  23.9  5.03e- 2 1.79e 0
#> 3 1       3         glm   3.08e 1     0.368     0.0342  31.7  2.18e- 2 1.57e 0
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3     25.5 
#>  2 1       1             20       14.1 
#>  3 1       1             10       12.3 
#>  4 1       1             23.3     10.9 
#>  5 1       1             10       41.7 
#>  6 1       1             10       39.1 
#>  7 1       1             20       17.8 
#>  8 1       1             10        8.46
#>  9 1       1             15       16.1 
#> 10 1       1             10        8.70
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
#> [1] "2026-10-21 13:39:14 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco + I(bio12^2) + I(elevation^2) + 
#>     I(sand^2) + bio12:elevation + bio12:sand + bio12:eco + elevation:sand + 
#>     eco:elevation + eco:sand
#> <environment: 0x6d427e3c0>
#> 
#> $metadata$sigma.formula
#> ~1
#> <environment: 0x6d427e3c0>
#> 
#> $metadata$nu.formula
#> ~1
#> <environment: 0x6d427e3c0>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6d427e3c0>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "ZAGA"
#> 
#> $metadata$hyperparameters$poly
#> [1] 2
#> 
#> $metadata$hyperparameters$inter_order
#> [1] 1
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_glm"
#> 
#> 

# Using third order polynomials and second order interaction terms
glm_3 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  distribution = "ZAGA",
  poly = 3,
  inter_order = 2,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco + I(bio12^2) + I(elevation^2) + I(sand^2) + I(bio12^3) + I(elevation^3) + I(sand^3) + bio12:elevation + bio12:sand + bio12:eco + elevation:sand + eco:elevation + eco:sand + bio12:elevation:sand + bio12:eco:elevation + bio12:eco:sand + eco:elevation:sand
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

glm_3
#> $model
#> 
#> Family:  c("ZAGA", "Zero Adjusted GA") 
#> Fitting method: RS() 
#> 
#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  
#>     nu.formula = nu_formula, tau.formula = tau_formula,  
#>     family = family, data = data, control = control_gamlss,      trace = FALSE) 
#> 
#> 
#> Mu Coefficients:
#>              (Intercept)                     bio12                 elevation  
#>                1.130e+03                -6.047e-01                 8.364e-01  
#>                     sand                     ecoAD                    ecoAMF  
#>               -4.644e-01                        NA                -1.125e+03  
#>                  ecoAPAF                   ecoBMDF                     ecoBY  
#>               -1.124e+03                        NA                        NA  
#>                     ecoC                   ecoCADP                    ecoCAP  
#>               -9.219e+02                        NA                -2.884e+03  
#>                   ecoCDF                     ecoCM                     ecoDC  
#>                       NA                        NA                -1.134e+03  
#>                     ecoE                     ecoHC                     ecoHM  
#>               -1.302e+02                -1.134e+03                        NA  
#>                    ecoHP                     ecoLM                      ecoP  
#>                       NA                        NA                        NA  
#>                   ecoPFS                     ecoPS                   ecoSABM  
#>               -1.226e+03                        NA                        NA  
#>                   ecoSAS                    ecoSAY                   ecoSCMS  
#>                       NA                -1.103e+03                        NA  
#>                 ecoSDMCF                     ecoUS                I(bio12^2)  
#>                       NA                        NA                -3.874e-06  
#>           I(elevation^2)                 I(sand^2)                I(bio12^3)  
#>                5.456e-06                 6.709e-06                 6.463e-10  
#>           I(elevation^3)                 I(sand^3)           bio12:elevation  
#>               -1.776e-09                -1.106e-09                 1.147e-05  
#>               bio12:sand               bio12:ecoAD              bio12:ecoAMF  
#>                3.008e-05                        NA                 6.110e-01  
#>            bio12:ecoAPAF             bio12:ecoBMDF               bio12:ecoBY  
#>                6.095e-01                        NA                        NA  
#>               bio12:ecoC             bio12:ecoCADP              bio12:ecoCAP  
#>                4.004e-01                        NA                        NA  
#>             bio12:ecoCDF               bio12:ecoCM               bio12:ecoDC  
#>                       NA                        NA                 6.137e-01  
#>               bio12:ecoE               bio12:ecoHC               bio12:ecoHM  
#>                       NA                 6.144e-01                        NA  
#>              bio12:ecoHP               bio12:ecoLM                bio12:ecoP  
#>                       NA                        NA                        NA  
#>             bio12:ecoPFS               bio12:ecoPS             bio12:ecoSABM  
#>                6.753e-01                        NA                        NA  
#>             bio12:ecoSAS              bio12:ecoSAY             bio12:ecoSCMS  
#>                       NA                 5.949e-01                        NA  
#>           bio12:ecoSDMCF               bio12:ecoUS            elevation:sand  
#>                       NA                        NA                 3.952e-05  
#>          elevation:ecoAD          elevation:ecoAMF         elevation:ecoAPAF  
#>                       NA                -8.626e-01                -8.548e-01  
#>        elevation:ecoBMDF           elevation:ecoBY            elevation:ecoC  
#>                       NA                        NA                -3.010e-02  
#>        elevation:ecoCADP          elevation:ecoCAP          elevation:ecoCDF  
#>                       NA                        NA                        NA  
#>          elevation:ecoCM           elevation:ecoDC            elevation:ecoE  
#>                       NA                -8.358e-01                        NA  
#>          elevation:ecoHC           elevation:ecoHM           elevation:ecoHP  
#>               -8.067e-01                        NA                        NA  
#>          elevation:ecoLM            elevation:ecoP          elevation:ecoPFS  
#>                       NA                        NA                -7.637e-01  
#>          elevation:ecoPS         elevation:ecoSABM          elevation:ecoSAS  
#>                       NA                        NA                        NA  
#>         elevation:ecoSAY         elevation:ecoSCMS        elevation:ecoSDMCF  
#>               -8.617e-01                        NA                        NA  
#>          elevation:ecoUS                sand:ecoAD               sand:ecoAMF  
#>                       NA                        NA                 4.550e-01  
#>             sand:ecoAPAF              sand:ecoBMDF                sand:ecoBY  
#>                4.447e-01                        NA                        NA  
#>                sand:ecoC              sand:ecoCADP               sand:ecoCAP  
#>                1.056e-01                        NA                        NA  
#>              sand:ecoCDF                sand:ecoCM                sand:ecoDC  
#>                       NA                        NA                 4.772e-01  
#>                sand:ecoE                sand:ecoHC                sand:ecoHM  
#>                       NA                 4.725e-01                        NA  
#>               sand:ecoHP                sand:ecoLM                 sand:ecoP  
#>                       NA                        NA                        NA  
#>              sand:ecoPFS                sand:ecoPS              sand:ecoSABM  
#>                4.057e-01                        NA                        NA  
#>              sand:ecoSAS               sand:ecoSAY              sand:ecoSCMS  
#>                       NA                 4.105e-01                        NA  
#>            sand:ecoSDMCF                sand:ecoUS      bio12:elevation:sand  
#>                       NA                        NA                -1.920e-08  
#>    bio12:elevation:ecoAD    bio12:elevation:ecoAMF   bio12:elevation:ecoAPAF  
#>                       NA                -3.754e-06                -5.410e-06  
#>  bio12:elevation:ecoBMDF     bio12:elevation:ecoBY      bio12:elevation:ecoC  
#>                       NA                        NA                -1.150e-04  
#>  bio12:elevation:ecoCADP    bio12:elevation:ecoCAP    bio12:elevation:ecoCDF  
#>                       NA                        NA                        NA  
#>    bio12:elevation:ecoCM     bio12:elevation:ecoDC      bio12:elevation:ecoE  
#>                       NA                -4.936e-06                        NA  
#>    bio12:elevation:ecoHC     bio12:elevation:ecoHM     bio12:elevation:ecoHP  
#>               -2.096e-05                        NA                        NA  
#>    bio12:elevation:ecoLM      bio12:elevation:ecoP    bio12:elevation:ecoPFS  
#>                       NA                        NA                        NA  
#>    bio12:elevation:ecoPS   bio12:elevation:ecoSABM    bio12:elevation:ecoSAS  
#>                       NA                        NA                        NA  
#>   bio12:elevation:ecoSAY   bio12:elevation:ecoSCMS  bio12:elevation:ecoSDMCF  
#>                       NA                        NA                        NA  
#>    bio12:elevation:ecoUS          bio12:sand:ecoAD         bio12:sand:ecoAMF  
#>                       NA                        NA                -3.330e-05  
#>       bio12:sand:ecoAPAF        bio12:sand:ecoBMDF          bio12:sand:ecoBY  
#>               -2.141e-05                        NA                        NA  
#>          bio12:sand:ecoC        bio12:sand:ecoCADP         bio12:sand:ecoCAP  
#>                3.161e-04                        NA                        NA  
#>        bio12:sand:ecoCDF          bio12:sand:ecoCM          bio12:sand:ecoDC  
#>                       NA                        NA                -3.649e-05  
#>          bio12:sand:ecoE          bio12:sand:ecoHC          bio12:sand:ecoHM  
#>                       NA                -3.658e-05                        NA  
#>         bio12:sand:ecoHP          bio12:sand:ecoLM           bio12:sand:ecoP  
#>                       NA                        NA                        NA  
#>        bio12:sand:ecoPFS          bio12:sand:ecoPS        bio12:sand:ecoSABM  
#>                       NA                        NA                        NA  
#>        bio12:sand:ecoSAS         bio12:sand:ecoSAY        bio12:sand:ecoSCMS  
#>                       NA                        NA                        NA  
#>      bio12:sand:ecoSDMCF          bio12:sand:ecoUS      elevation:sand:ecoAD  
#>                       NA                        NA                        NA  
#>    elevation:sand:ecoAMF    elevation:sand:ecoAPAF    elevation:sand:ecoBMDF  
#>                2.533e-05                -6.096e-06                        NA  
#>     elevation:sand:ecoBY       elevation:sand:ecoC    elevation:sand:ecoCADP  
#>                       NA                -9.980e-04                        NA  
#>    elevation:sand:ecoCAP     elevation:sand:ecoCDF      elevation:sand:ecoCM  
#>                       NA                        NA                        NA  
#>     elevation:sand:ecoDC       elevation:sand:ecoE      elevation:sand:ecoHC  
#>               -5.465e-05                        NA                -6.337e-05  
#>     elevation:sand:ecoHM      elevation:sand:ecoHP      elevation:sand:ecoLM  
#>                       NA                        NA                        NA  
#>      elevation:sand:ecoP     elevation:sand:ecoPFS      elevation:sand:ecoPS  
#>                       NA                        NA                        NA  
#>   elevation:sand:ecoSABM     elevation:sand:ecoSAS     elevation:sand:ecoSAY  
#>                       NA                        NA                        NA  
#>   elevation:sand:ecoSCMS   elevation:sand:ecoSDMCF      elevation:sand:ecoUS  
#>                       NA                        NA                        NA  
#> Sigma Coefficients:
#> (Intercept)  
#>     -0.2717  
#> Nu Coefficients:
#> (Intercept)  
#>      -1.605  
#> 
#>  Degrees of Freedom for the fit: 61 Residual Deg. of Freedom   501 
#> Global Deviance:     4520.12 
#>             AIC:     4642.12 
#>             SBC:     4906.34 
#> 
#> $predictors
#> # A tibble: 1 × 6
#>   model response c1    c2        c3    f    
#>   <chr> <chr>    <chr> <chr>     <chr> <chr>
#> 1 glm   ind_ha   bio12 elevation sand  eco  
#> 
#> $performance
#> # A tibble: 1 × 13
#>   model mae_mean  mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>   <chr>    <dbl>   <dbl>           <dbl>         <dbl>          <dbl>
#> 1 glm    1.56e29 2.70e29           0.221         0.149        -0.0191
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model     mae corr_spear corr_pear inter     slope   pdisp
#>   <chr>   <chr>     <chr>   <dbl>      <dbl>     <dbl> <dbl>     <dbl>   <dbl>
#> 1 1       1         glm   1.70e24     0.0539   -0.0579  27.1 -8.29e-26 6.99e23
#> 2 1       2         glm   3.37e 1     0.268     0.0233  25.4  7.13e- 3 3.26e 0
#> 3 1       3         glm   4.68e29     0.341    -0.0226  32.6 -1.53e-31 1.48e29
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3     24.0 
#>  2 1       1             20       19.5 
#>  3 1       1             10       13.4 
#>  4 1       1             23.3      7.53
#>  5 1       1             10       44.5 
#>  6 1       1             10       43.8 
#>  7 1       1             20       16.7 
#>  8 1       1             10        8.97
#>  9 1       1             15       16.2 
#> 10 1       1             10        8.67
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
#> [1] "2026-10-21 13:39:14 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco + I(bio12^2) + I(elevation^2) + 
#>     I(sand^2) + I(bio12^3) + I(elevation^3) + I(sand^3) + bio12:elevation + 
#>     bio12:sand + bio12:eco + elevation:sand + eco:elevation + 
#>     eco:sand + bio12:elevation:sand + bio12:eco:elevation + bio12:eco:sand + 
#>     eco:elevation:sand
#> <environment: 0x6df1cc2a8>
#> 
#> $metadata$sigma.formula
#> ~1
#> <environment: 0x6df1cc2a8>
#> 
#> $metadata$nu.formula
#> ~1
#> <environment: 0x6df1cc2a8>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6df1cc2a8>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "ZAGA"
#> 
#> $metadata$hyperparameters$poly
#> [1] 3
#> 
#> $metadata$hyperparameters$inter_order
#> [1] 2
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_glm"
#> 
#> 

# Setting formulas for different distribution parameters
glm_4 <- fit_abund_glm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("bio12", "elevation", "sand"),
  predictors_f = c("eco"),
  partition = ".part1",
  distribution = "ZAGA",
  fit_formula = ind_ha ~ bio12 + elevation + sand + eco,
  sigma_formula = ind_ha ~ bio12 + elevation + sand,
  poly = 0,
  inter_order = 0,
  predict_part = TRUE
)
#> Formula used for model fitting:
#> ind_ha ~ bio12 + elevation + sand + eco
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

glm_4
#> $model
#> 
#> Family:  c("ZAGA", "Zero Adjusted GA") 
#> Fitting method: RS() 
#> 
#> Call:  gamlss::gamlss(formula = formula1, sigma.formula = sigma_formula,  
#>     nu.formula = nu_formula, tau.formula = tau_formula,  
#>     family = family, data = data, control = control_gamlss,      trace = FALSE) 
#> 
#> 
#> Mu Coefficients:
#> (Intercept)        bio12    elevation         sand        ecoAD       ecoAMF  
#>   5.5976032   -0.0004723   -0.0008510   -0.0003508           NA   -1.3940197  
#>     ecoAPAF      ecoBMDF        ecoBY         ecoC      ecoCADP       ecoCAP  
#>  -1.4039557           NA           NA   -1.2443525           NA    1.6689538  
#>      ecoCDF        ecoCM        ecoDC         ecoE        ecoHC        ecoHM  
#>          NA           NA   -1.3775897   -2.0709633   -0.8742266           NA  
#>       ecoHP        ecoLM         ecoP       ecoPFS        ecoPS      ecoSABM  
#>          NA           NA           NA   -1.8515159           NA           NA  
#>      ecoSAS       ecoSAY      ecoSCMS     ecoSDMCF        ecoUS  
#>          NA   -0.5340199           NA           NA           NA  
#> Sigma Coefficients:
#> (Intercept)        bio12    elevation         sand  
#>  -2.061e-01   -9.656e-05   -4.136e-05    5.685e-04  
#> Nu Coefficients:
#> (Intercept)  
#>      -1.605  
#> 
#>  Degrees of Freedom for the fit: 18 Residual Deg. of Freedom   544 
#> Global Deviance:     4578.26 
#>             AIC:     4614.26 
#>             SBC:     4692.23 
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
#> 1 glm       25.8   1.60           0.305         0.121          0.279
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
#> 
#> $performance_part
#> # A tibble: 3 × 9
#>   replica partition model   mae corr_spear corr_pear inter slope pdisp
#>   <chr>   <chr>     <chr> <dbl>      <dbl>     <dbl> <dbl> <dbl> <dbl>
#> 1 1       1         glm    27.3      0.222     0.188 15.3  0.313 0.601
#> 2 1       2         glm    24.1      0.250     0.242  8.69 0.517 0.468
#> 3 1       3         glm    25.9      0.444     0.406 -8.29 1.27  0.321
#> 
#> $predicted_part
#> # A tibble: 562 × 4
#>    replica partition observed predicted
#>    <chr>   <chr>        <dbl>     <dbl>
#>  1 1       1             13.3      37.3
#>  2 1       1             20        15.5
#>  3 1       1             10        11.6
#>  4 1       1             23.3      12.3
#>  5 1       1             10        28.8
#>  6 1       1             10        26.3
#>  7 1       1             20        19.8
#>  8 1       1             10        15.1
#>  9 1       1             15        17.4
#> 10 1       1             10        13.5
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
#> [1] "2026-10-21 13:39:14 PDT"
#> 
#> $metadata$formula
#> ind_ha ~ bio12 + elevation + sand + eco
#> <environment: 0x6f4948a50>
#> 
#> $metadata$sigma.formula
#> ind_ha ~ bio12 + elevation + sand
#> <environment: 0x6f4948a50>
#> 
#> $metadata$nu.formula
#> ~1
#> <environment: 0x6df7ca158>
#> 
#> $metadata$tau.formula
#> ~1
#> <environment: 0x6df7ca158>
#> 
#> $metadata$hyperparameters
#> $metadata$hyperparameters$distribution
#> [1] "ZAGA"
#> 
#> $metadata$hyperparameters$poly
#> [1] 0
#> 
#> $metadata$hyperparameters$inter_order
#> [1] 0
#> 
#> 
#> $metadata$source_function
#> [1] "fit_abund_glm"
#> 
#> 
# }
```

# Merge model performance tables

Merge model performance tables

## Usage

``` r
adm_summarize(models)
```

## Arguments

- models:

  list. A list of a single or several models fitted with some of fit\_
  or tune\_ functions. Usage models = list(mod1, mod2, mod3)

## Value

A tibble object with combined model performance for all input models.
Models fit with tune will include model performance for the best
hyperparameters.

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)

data("sppabund")
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Species abundance data, coordinates, and partition
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y, .part1)
some_sp
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y .part1
#>    <chr>        <dbl> <dbl> <dbl>  <int>
#>  1 Species one   10   -61.2 -25.0      2
#>  2 Species one   10   -60.7 -25.0      3
#>  3 Species one   10   -60.2 -25.0      2
#>  4 Species one   13.3 -60.7 -25.4      1
#>  5 Species one   10   -60.2 -25.9      3
#>  6 Species one   10   -54.0 -26.0      3
#>  7 Species one   20   -53.9 -26.0      1
#>  8 Species one   10   -53.9 -26.0      1
#>  9 Species one   10   -53.9 -26.1      3
#> 10 Species one   20   -53.9 -26.1      2
#> # ℹ 928 more rows

# Extract data
some_sp <-
  adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar
  )

# Fit RAF
m_raf <- fit_abund_raf(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part1",
)
#> Formula used for model fitting:
#> ind_ha ~ elevation + sand + bio3 + bio12
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

# Fit SVM
m_svm <- fit_abund_svm(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part1"
)
#> Formula used for model fitting:
#> ind_ha ~ elevation + sand + bio3 + bio12
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3

# XGB
m_xbg <- fit_abund_xgb(
  data = some_sp,
  response = "ind_ha",
  predictors = c("elevation", "sand", "bio3", "bio12"),
  partition = ".part1"
)
#> Replica number: 1/1
#> -- Partition number 1/3
#> -- Partition number 2/3
#> -- Partition number 3/3


perf <- adm_summarize(list(m_svm, m_raf, m_xbg))

perf
#> # A tibble: 3 × 14
#>   model_ID model mae_mean mae_sd corr_spear_mean corr_spear_sd corr_pear_mean
#>      <int> <chr>    <dbl>  <dbl>           <dbl>         <dbl>          <dbl>
#> 1        1 svm       14.4  0.793           0.616        0.0774          0.450
#> 2        2 raf       15.9  0.562           0.594        0.0632          0.417
#> 3        3 xgb       16.7  0.177           0.581        0.0829          0.397
#> # ℹ 7 more variables: corr_pear_sd <dbl>, inter_mean <dbl>, inter_sd <dbl>,
#> #   slope_mean <dbl>, slope_sd <dbl>, pdisp_mean <dbl>, pdisp_sd <dbl>
# }
```

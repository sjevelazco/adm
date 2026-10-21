# Performs data transformation on a variable based on the specified method.

This function transforms the data in a tibble or SpatRaster object on
the method specified The available methods are "01", "zscore", "log",
and "round".

## Usage

``` r
adm_transform(data, variable, method, inverse = FALSE, t_terms = NULL)
```

## Arguments

- data:

  A data.frame, tibble, or SpatRaster containing the data.

- variable:

  A character string specifying the variable (column) to be transformed.

- method:

  A character string specifying the method to be used for
  transformation. Available methods are "01", "zscore", "log", and
  "round."

  - For "01", it scales the variable between 0 and 1 using the formula
    (x - min(x)) (max(x) - min(x)).

  - For "zscore", it standardizes the variable by subtracting the mean
    and dividing by the standard deviation.

  - For "log", it applies natural logarithm transformation to the
    variable.

  - For "log1", it sums 1 and then applies natural logarithm
    transformation to the variable.

  - For "round", it rounds the variable's values to the nearest whole
    numbers.

- inverse:

  logical. Invert the transformation?

- t_terms:

  vector. c(a,b):

  - For "01", a = min(x), b = max(x).

  - For "zscore", a = mean(x), b = sd(x).

  - For "log" and "log1, not needed.

  - Can't invert "round" transformations.

## Value

A data.frame or tibble with the transformed variable added as a new
column. The new column's name is the original variable name followed by
an underscore and method name.

## Examples

``` r
# \donttest{
require(dplyr)
require(terra)

# Select data for a single species
data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y)

envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)[["bio12"]]

# Transform tabular data
## Transform abundance data to 0-1
some_sp_2 <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "01"
)
some_sp_2
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_01
#>    <chr>        <dbl> <dbl> <dbl>     <dbl>
#>  1 Species one   10   -61.2 -25.0    0.0351
#>  2 Species one   10   -60.7 -25.0    0.0351
#>  3 Species one   10   -60.2 -25.0    0.0351
#>  4 Species one   13.3 -60.7 -25.4    0.0467
#>  5 Species one   10   -60.2 -25.9    0.0351
#>  6 Species one   10   -54.0 -26.0    0.0351
#>  7 Species one   20   -53.9 -26.0    0.0701
#>  8 Species one   10   -53.9 -26.0    0.0351
#>  9 Species one   10   -53.9 -26.1    0.0351
#> 10 Species one   20   -53.9 -26.1    0.0701
#> # ℹ 928 more rows

## Transform abundance data z-score
some_sp_2 <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "zscore"
)
some_sp_2
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_zscore
#>    <chr>        <dbl> <dbl> <dbl>         <dbl>
#>  1 Species one   10   -61.2 -25.0       -0.218 
#>  2 Species one   10   -60.7 -25.0       -0.218 
#>  3 Species one   10   -60.2 -25.0       -0.218 
#>  4 Species one   13.3 -60.7 -25.4       -0.113 
#>  5 Species one   10   -60.2 -25.9       -0.218 
#>  6 Species one   10   -54.0 -26.0       -0.218 
#>  7 Species one   20   -53.9 -26.0        0.0961
#>  8 Species one   10   -53.9 -26.0       -0.218 
#>  9 Species one   10   -53.9 -26.1       -0.218 
#> 10 Species one   20   -53.9 -26.1        0.0961
#> # ℹ 928 more rows

## Transform abundance data log
some_sp_2 <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "log"
)
some_sp_2
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_log
#>    <chr>        <dbl> <dbl> <dbl>      <dbl>
#>  1 Species one   10   -61.2 -25.0       2.30
#>  2 Species one   10   -60.7 -25.0       2.30
#>  3 Species one   10   -60.2 -25.0       2.30
#>  4 Species one   13.3 -60.7 -25.4       2.59
#>  5 Species one   10   -60.2 -25.9       2.30
#>  6 Species one   10   -54.0 -26.0       2.30
#>  7 Species one   20   -53.9 -26.0       3.00
#>  8 Species one   10   -53.9 -26.0       2.30
#>  9 Species one   10   -53.9 -26.1       2.30
#> 10 Species one   20   -53.9 -26.1       3.00
#> # ℹ 928 more rows

## Round abundance data
some_sp_2 <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "round"
)
some_sp_2
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_round
#>    <chr>        <dbl> <dbl> <dbl>        <dbl>
#>  1 Species one   10   -61.2 -25.0           10
#>  2 Species one   10   -60.7 -25.0           10
#>  3 Species one   10   -60.2 -25.0           10
#>  4 Species one   13.3 -60.7 -25.4           13
#>  5 Species one   10   -60.2 -25.9           10
#>  6 Species one   10   -54.0 -26.0           10
#>  7 Species one   20   -53.9 -26.0           20
#>  8 Species one   10   -53.9 -26.0           10
#>  9 Species one   10   -53.9 -26.1           10
#> 10 Species one   20   -53.9 -26.1           20
#> # ℹ 928 more rows

# Tranform raster data
## Transform to 0-1
envar_2 <- adm_transform(
  data = envar,
  variable = "bio12",
  method = "01"
)
envar_2
#> class       : SpatRaster
#> size        : 201, 274, 2  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory
#> varnames    : envar
#>               envar
#> names       :       bio12, bio12_01
#> min values  :       9.283,        0
#> max values  : 3171.336914,        1

## Transform z-score
envar_2 <- adm_transform(
  data = envar,
  variable = "bio12",
  method = "zscore"
)
envar_2
#> class       : SpatRaster
#> size        : 201, 274, 2  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory
#> varnames    : envar
#>               envar
#> names       :       bio12, bio12_zscore
#> min values  :       9.283,    -2.249803
#> max values  : 3171.336914,     3.899979

## Transform log
envar_2 <- adm_transform(
  data = envar,
  variable = "bio12",
  method = "log"
)
envar_2
#> class       : SpatRaster
#> size        : 201, 274, 2  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory
#> varnames    : envar
#>               envar
#> names       :       bio12, bio12_log
#> min values  :       9.283,  2.228185
#> max values  : 3171.336914,  8.061909

## Round
envar_2 <- adm_transform(
  data = envar,
  variable = "bio12",
  method = "round"
)
envar_2
#> class       : SpatRaster
#> size        : 201, 274, 2  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory
#> varnames    : envar
#>               envar
#> names       :       bio12, bio12_round
#> min values  :       9.283,           9
#> max values  : 3171.336914,        3171

# Invert transformation
## Invert 01 tranformation
some_sp_transformed <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "01"
)
some_sp_transformed
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_01
#>    <chr>        <dbl> <dbl> <dbl>     <dbl>
#>  1 Species one   10   -61.2 -25.0    0.0351
#>  2 Species one   10   -60.7 -25.0    0.0351
#>  3 Species one   10   -60.2 -25.0    0.0351
#>  4 Species one   13.3 -60.7 -25.4    0.0467
#>  5 Species one   10   -60.2 -25.9    0.0351
#>  6 Species one   10   -54.0 -26.0    0.0351
#>  7 Species one   20   -53.9 -26.0    0.0701
#>  8 Species one   10   -53.9 -26.0    0.0351
#>  9 Species one   10   -53.9 -26.1    0.0351
#> 10 Species one   20   -53.9 -26.1    0.0701
#> # ℹ 928 more rows

some_sp_inverted <- adm_transform(
  data = some_sp_transformed,
  variable = "ind_ha_01",
  method = "01",
  inverse = TRUE,
  t_terms = c(
    a = min(some_sp[["ind_ha"]]),
    b = max(some_sp[["ind_ha"]])
  )
)
some_sp_inverted
#> # A tibble: 938 × 6
#>    species     ind_ha     x     y ind_ha_01 ind_ha_01_inverted
#>    <chr>        <dbl> <dbl> <dbl>     <dbl>              <dbl>
#>  1 Species one   10   -61.2 -25.0    0.0351               10  
#>  2 Species one   10   -60.7 -25.0    0.0351               10  
#>  3 Species one   10   -60.2 -25.0    0.0351               10  
#>  4 Species one   13.3 -60.7 -25.4    0.0467               13.3
#>  5 Species one   10   -60.2 -25.9    0.0351               10  
#>  6 Species one   10   -54.0 -26.0    0.0351               10  
#>  7 Species one   20   -53.9 -26.0    0.0701               20  
#>  8 Species one   10   -53.9 -26.0    0.0351               10  
#>  9 Species one   10   -53.9 -26.1    0.0351               10  
#> 10 Species one   20   -53.9 -26.1    0.0701               20  
#> # ℹ 928 more rows

## Invert z-score tranformation
some_sp_transformed <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "zscore"
)
some_sp_transformed
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_zscore
#>    <chr>        <dbl> <dbl> <dbl>         <dbl>
#>  1 Species one   10   -61.2 -25.0       -0.218 
#>  2 Species one   10   -60.7 -25.0       -0.218 
#>  3 Species one   10   -60.2 -25.0       -0.218 
#>  4 Species one   13.3 -60.7 -25.4       -0.113 
#>  5 Species one   10   -60.2 -25.9       -0.218 
#>  6 Species one   10   -54.0 -26.0       -0.218 
#>  7 Species one   20   -53.9 -26.0        0.0961
#>  8 Species one   10   -53.9 -26.0       -0.218 
#>  9 Species one   10   -53.9 -26.1       -0.218 
#> 10 Species one   20   -53.9 -26.1        0.0961
#> # ℹ 928 more rows

some_sp_inverted <- adm_transform(
  data = some_sp_transformed,
  variable = "ind_ha_zscore",
  method = "zscore",
  inverse = TRUE,
  t_terms = c(
    a = mean(some_sp[["ind_ha"]]),
    b = sd(some_sp[["ind_ha"]])
  )
)
some_sp_inverted
#> # A tibble: 938 × 6
#>    species     ind_ha     x     y ind_ha_zscore ind_ha_zscore_inverted
#>    <chr>        <dbl> <dbl> <dbl>         <dbl>                  <dbl>
#>  1 Species one   10   -61.2 -25.0       -0.218                    10  
#>  2 Species one   10   -60.7 -25.0       -0.218                    10  
#>  3 Species one   10   -60.2 -25.0       -0.218                    10  
#>  4 Species one   13.3 -60.7 -25.4       -0.113                    13.3
#>  5 Species one   10   -60.2 -25.9       -0.218                    10  
#>  6 Species one   10   -54.0 -26.0       -0.218                    10  
#>  7 Species one   20   -53.9 -26.0        0.0961                   20  
#>  8 Species one   10   -53.9 -26.0       -0.218                    10  
#>  9 Species one   10   -53.9 -26.1       -0.218                    10  
#> 10 Species one   20   -53.9 -26.1        0.0961                   20  
#> # ℹ 928 more rows

## Invert log and log1
some_sp_transformed <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "log"
)
some_sp_transformed
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_log
#>    <chr>        <dbl> <dbl> <dbl>      <dbl>
#>  1 Species one   10   -61.2 -25.0       2.30
#>  2 Species one   10   -60.7 -25.0       2.30
#>  3 Species one   10   -60.2 -25.0       2.30
#>  4 Species one   13.3 -60.7 -25.4       2.59
#>  5 Species one   10   -60.2 -25.9       2.30
#>  6 Species one   10   -54.0 -26.0       2.30
#>  7 Species one   20   -53.9 -26.0       3.00
#>  8 Species one   10   -53.9 -26.0       2.30
#>  9 Species one   10   -53.9 -26.1       2.30
#> 10 Species one   20   -53.9 -26.1       3.00
#> # ℹ 928 more rows

some_sp_inverted <- adm_transform(
  data = some_sp_transformed,
  variable = "ind_ha_log",
  method = "log",
  inverse = TRUE
)
#> Warning: log and log1 tranformations can't restore original negative values, use carefully.
some_sp_inverted
#> # A tibble: 938 × 6
#>    species     ind_ha     x     y ind_ha_log ind_ha_log_inverted
#>    <chr>        <dbl> <dbl> <dbl>      <dbl>               <dbl>
#>  1 Species one   10   -61.2 -25.0       2.30                10  
#>  2 Species one   10   -60.7 -25.0       2.30                10  
#>  3 Species one   10   -60.2 -25.0       2.30                10  
#>  4 Species one   13.3 -60.7 -25.4       2.59                13.3
#>  5 Species one   10   -60.2 -25.9       2.30                10  
#>  6 Species one   10   -54.0 -26.0       2.30                10  
#>  7 Species one   20   -53.9 -26.0       3.00                20  
#>  8 Species one   10   -53.9 -26.0       2.30                10  
#>  9 Species one   10   -53.9 -26.1       2.30                10  
#> 10 Species one   20   -53.9 -26.1       3.00                20  
#> # ℹ 928 more rows

some_sp_transformed <- adm_transform(
  data = some_sp,
  variable = "ind_ha",
  method = "log1"
)
some_sp_transformed
#> # A tibble: 938 × 5
#>    species     ind_ha     x     y ind_ha_log1
#>    <chr>        <dbl> <dbl> <dbl>       <dbl>
#>  1 Species one   10   -61.2 -25.0        2.40
#>  2 Species one   10   -60.7 -25.0        2.40
#>  3 Species one   10   -60.2 -25.0        2.40
#>  4 Species one   13.3 -60.7 -25.4        2.66
#>  5 Species one   10   -60.2 -25.9        2.40
#>  6 Species one   10   -54.0 -26.0        2.40
#>  7 Species one   20   -53.9 -26.0        3.04
#>  8 Species one   10   -53.9 -26.0        2.40
#>  9 Species one   10   -53.9 -26.1        2.40
#> 10 Species one   20   -53.9 -26.1        3.04
#> # ℹ 928 more rows

some_sp_inverted <- adm_transform(
  data = some_sp_transformed,
  variable = "ind_ha_log1",
  method = "log1",
  inverse = TRUE
)
#> Warning: log and log1 tranformations can't restore original negative values, use carefully.
some_sp_inverted
#> # A tibble: 938 × 6
#>    species     ind_ha     x     y ind_ha_log1 ind_ha_log1_inverted
#>    <chr>        <dbl> <dbl> <dbl>       <dbl>                <dbl>
#>  1 Species one   10   -61.2 -25.0        2.40                 10  
#>  2 Species one   10   -60.7 -25.0        2.40                 10  
#>  3 Species one   10   -60.2 -25.0        2.40                 10  
#>  4 Species one   13.3 -60.7 -25.4        2.66                 13.3
#>  5 Species one   10   -60.2 -25.9        2.40                 10  
#>  6 Species one   10   -54.0 -26.0        2.40                 10  
#>  7 Species one   20   -53.9 -26.0        3.04                 20  
#>  8 Species one   10   -53.9 -26.0        2.40                 10  
#>  9 Species one   10   -53.9 -26.1        2.40                 10  
#> 10 Species one   20   -53.9 -26.1        3.04                 20  
#> # ℹ 928 more rows

## To invert raster, is the same process
## Example:
envar_transformed <- adm_transform(
  data = envar,
  variable = "bio12",
  method = "01"
)
envar_transformed
#> class       : SpatRaster
#> size        : 201, 274, 2  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory
#> varnames    : envar
#>               envar
#> names       :       bio12, bio12_01
#> min values  :       9.283,        0
#> max values  : 3171.336914,        1

envar_inverted <- adm_transform(
  data = envar_transformed,
  variable = "bio12_01",
  method = "01",
  inverse = TRUE,
  t_terms = c(
    a = terra::global(envar, min, na.rm = TRUE),
    b = terra::global(envar, max, na.rm = TRUE)
  )
)
envar_inverted
#> class       : SpatRaster
#> size        : 201, 274, 3  (nrow, ncol, nlyr)
#> resolution  : 0.08333333, 0.08333333  (x, y)
#> extent      : -69.41681, -46.58347, -34.38347, -17.63347  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326)
#> sources     : envar.tif
#>               memory (2 layers)
#> varnames    : envar
#>               envar
#>               envar
#> names       :       bio12, bio12_01, bio12_01_inverted
#> min values  :       9.283,        0,             9.283
#> max values  : 3171.336914,        1,       3171.336914
# }
```

# Extract values from a spatial raster based on x and y coordinates

The function extracts environmental data at the given x and y
coordinates

## Usage

``` r
adm_extract(data, x, y, env_layer, variables = NULL, filter_na = TRUE)
```

## Arguments

- data:

  data.frame or tibble. Database with species abundance with x and y
  coordinates

- x:

  character. Column name with spatial x coordinates

- y:

  character. Column name with spatial y coordinates

- env_layer:

  SpatRaster. Raster with environmental variables.

- variables:

  character. Vector with the variable names of predictor (environmental)
  variables Usage variables = c("elevation", "sand", "cfvo"). If no
  variable is specified, function will return data for all layers.
  Default NULL

- filter_na:

  logical. If filter_na = TRUE (default), the rows with NA values for
  any of the environmental variables are removed from the returned
  tibble.

## Value

A tibble that returns the original data base with additional columns for
the extracted environmental variables at each x y location from the
SpatRaster object used in 'env_layer'

## Examples

``` r
# \donttest{
require(dplyr)
#> Loading required package: dplyr
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
require(terra)
#> Loading required package: terra
#> terra 1.9.50
#> 
#> Attaching package: ‘terra’
#> The following objects are masked from ‘package:testthat’:
#> 
#>     compare, describe

# Datasbase with species abundance and x and y coordinates
data("sppabund")

# Raster data with environmental variables
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)

# Extract data for a single species
some_sp <- sppabund %>%
  filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y)

# Extract environmental data from envar raster for all locations in spp
ex_spp <-
  adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar,
    variables = NULL,
    filter_na = FALSE
  )

# Extract environmental for two variables and remove rows with NAs
ex_spp2 <-
  adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar,
    variables = c("bio1", "elevation"),
    filter_na = TRUE
  )

ex_spp
#> # A tibble: 938 × 12
#>    species     ind_ha     x     y  bio1 bio12 bio15  bio3  cfvo elevation  sand
#>    <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>     <dbl> <dbl>
#>  1 Species one   10   -61.2 -25.0  22.9  921.  66.7 0.475  33.5      144.  217.
#>  2 Species one   10   -60.7 -25.0  22.8 1008.  62.5 0.461  31.9      134.  163.
#>  3 Species one   10   -60.2 -25.0  22.7 1103.  55.8 0.458  26.2      118.  185.
#>  4 Species one   13.3 -60.7 -25.4  22.5 1056.  61.0 0.457  41.4      124.  162.
#>  5 Species one   10   -60.2 -25.9  22.1 1175.  55.8 0.449  43.6      104.  213.
#>  6 Species one   10   -54.0 -26.0  20.7 2185.  19.0 0.494  60.2      411.  125.
#>  7 Species one   20   -53.9 -26.0  20.7 2185.  19.0 0.494  60.2      411.  125.
#>  8 Species one   10   -53.9 -26.0  20.3 2403.  18.7 0.497  53.4      488.  113.
#>  9 Species one   10   -53.9 -26.1  20.0 2349.  18.3 0.496  55.5      509.  117.
#> 10 Species one   20   -53.9 -26.1  20.0 2349.  18.3 0.496  55.5      509.  117.
#> # ℹ 928 more rows
#> # ℹ 1 more variable: eco <fct>
ex_spp2
#> # A tibble: 938 × 6
#>    species     ind_ha     x     y  bio1 elevation
#>    <chr>        <dbl> <dbl> <dbl> <dbl>     <dbl>
#>  1 Species one   10   -61.2 -25.0  22.9      144.
#>  2 Species one   10   -60.7 -25.0  22.8      134.
#>  3 Species one   10   -60.2 -25.0  22.7      118.
#>  4 Species one   13.3 -60.7 -25.4  22.5      124.
#>  5 Species one   10   -60.2 -25.9  22.1      104.
#>  6 Species one   10   -54.0 -26.0  20.7      411.
#>  7 Species one   20   -53.9 -26.0  20.7      411.
#>  8 Species one   10   -53.9 -26.0  20.3      488.
#>  9 Species one   10   -53.9 -26.1  20.0      509.
#> 10 Species one   20   -53.9 -26.1  20.0      509.
#> # ℹ 928 more rows
# }
```

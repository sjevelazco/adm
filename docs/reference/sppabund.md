# A data set containing species abundance of three species, partition folds, and environmental variables.

A data set containing species abundance of three species, partition
folds, and environmental variables.

## Usage

``` r
sppabund
```

## Format

A tibble with 2767 rows and 12 variables:

- species:

  species names

- ind_ha:

  species abundance expressed as individuals per hectare

- x:

  longitude of species occurrences

- y:

  latitude of species occurrences

- bio1:

  bioclimatic variable related to annual mean temperature

- bio3:

  bioclimatic variable related to isothermality

- bio12:

  bioclimatic variable related to annual precipitation

- bio15:

  bioclimatic variable related to precipitation seasonality

- cfvo:

  edaphic variable related to volumetric fraction of coarse fragments

- elevation:

  topographic elevation

- sand:

  edaphic variable related to soil sand content

- eoc:

  ecoregion

- .part1 ... .part3:

  repeate k-folds

## Examples

``` r
# \donttest{
require(dplyr)
data("sppabund")
sppabund
#> # A tibble: 2,767 × 16
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
#> # ℹ 2,757 more rows
#> # ℹ 5 more variables: eco <fct>, pr_ab <dbl>, .part1 <int>, .part2 <int>,
#> #   .part3 <int>
# }
```

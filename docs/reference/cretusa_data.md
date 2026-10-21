# A data set containing abundance of Cynophalla retusa.

Cynophalla retusa (Griseb.) Cornejo & Iltis (Capparaceae) is a shrub
native to northeastern Argentina, Paraguay, Bolivia, and central Brazil.

## Usage

``` r
cretusa_data
```

## Format

A tibble with 366 rows and 6 variables:

- species:

  species names

- ind_ha:

  species abundance expressed as individuals per hectare

- x:

  longitude of species occurrences

- y:

  latitude of species occurrences

- .part:

  partition folds

## Examples

``` r
# \donttest{
require(dplyr)
data("cretusa_data")
cretusa_data
#> # A tibble: 366 × 5
#>    species           ind_ha     x     y .part
#>    <chr>              <int> <dbl> <dbl> <int>
#>  1 Cynophalla retusa     10 -64.5 -22.7     1
#>  2 Cynophalla retusa     10 -64.1 -22.7     1
#>  3 Cynophalla retusa     20 -64.7 -23.1     2
#>  4 Cynophalla retusa      0 -62.6 -23.2     3
#>  5 Cynophalla retusa      0 -61.7 -24.5     3
#>  6 Cynophalla retusa      0 -61.6 -25.0     3
#>  7 Cynophalla retusa      0 -61.2 -24.5     3
#>  8 Cynophalla retusa      0 -64.9 -23.8     2
#>  9 Cynophalla retusa      0 -65.3 -24.4     3
#> 10 Cynophalla retusa      0 -64.7 -24.8     1
#> # ℹ 356 more rows
# }
```

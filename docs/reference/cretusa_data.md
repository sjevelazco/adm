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
if (FALSE) { # \dontrun{
require(dplyr)
data("cretusa_data")
cretusa_data
} # }
```

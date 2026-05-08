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
if (FALSE) { # \dontrun{
require(dplyr)
data("sppabund")
sppabund
} # }
```

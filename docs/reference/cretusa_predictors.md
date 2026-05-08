# Raster with Principal Component

A raster with principal components derived by a principal component
analysis based on climatic (from Chelsa: chelsa-climate.org) and edaphic
(from SoilGrids: soilgrids.org) variables.

## Usage

``` r
cretusa_predictors()
```

## Format

A raster in tif format with the first principal components.

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
envar <- system.file("external/cretusa_predictors.tif", package = "adm")
envar <- terra::rast(envar)
plot(envar)
} # }
```

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

## Value

A `SpatRaster` (see
[`rast`](https://rspatial.github.io/terra/reference/rast.html)) with the
first principal components (PC1, PC2, ...) read from the raster file
shipped with the package (`external/cretusa_predictors.tif`).

## Examples

``` r
# \donttest{
require(terra)
predictors <- adm::cretusa_predictors()
plot(predictors)

# }
```

# Raster with environmental data

A raster with climatic (from Chelsa: chelsa-climate.org), edaphic (from
SoilGrids: soilgrids.org), and ecoregion data (from WWF
https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world).

## Usage

``` r
envar()
```

## Format

A raster in tif format with the following variables.

- bio1:

  Annual mean temperature

- bio12:

  Annual precipitation

- bio15:

  Precipitation seasonality

- bio3:

  Isothermality

- cfvo:

  Volumetric fraction of coarse fragments

- elevation:

  Elevation

- sand:

  Soil sand content

- eco:

  Terrestrial ecoregion

## Examples

``` r
if (FALSE) { # \dontrun{
require(terra)
envar <- system.file("external/envar.tif", package = "adm")
envar <- terra::rast(envar)
plot(envar)
} # }
```

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

## Value

A `SpatRaster` (see
[`rast`](https://rspatial.github.io/terra/reference/rast.html)) with
eight layers (bio1, bio12, bio15, bio3, cfvo, elevation, sand, and eco)
read from the raster file shipped with the package
(`external/envar.tif`).

## Examples

``` r
# \donttest{
require(terra)
envar <- adm::envar()
plot(envar)

# }
```

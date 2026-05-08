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
if (FALSE) { # \dontrun{
require(terra)

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
ex_spp2
} # }
```

# Small shared fixtures for the lightweight tests (kept tiny to be CRAN-friendly)
suppressMessages({
  library(dplyr)
  library(terra)
})

data("sppabund", package = "adm")
test_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)
test_sp <- balance_dataset(test_sp, response = "ind_ha", absence_ratio = 0.2)
test_pred <- c("bio12", "elevation", "sand")
test_envar <- terra::rast(system.file("external/envar.tif", package = "adm"))

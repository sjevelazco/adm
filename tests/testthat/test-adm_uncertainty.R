test_that("adm_uncertainty returns a raster", {
  m <- fit_abund_raf(test_sp, "ind_ha", test_pred,
    partition = ".part1", ntree = 20, verbose = FALSE
  )
  u <- suppressWarnings(adm_uncertainty(m, test_sp, "ind_ha", test_envar, iteration = 2))
  expect_s4_class(u, "SpatRaster")
  expect_equal(terra::ncell(u), terra::ncell(test_envar))
  expect_true(all(terra::global(u, "min", na.rm = TRUE)[, 1] >= 0))
})

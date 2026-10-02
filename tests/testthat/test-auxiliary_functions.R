test_that("res_calculate", {
  expect_equal(res_calculate("layer", in_res = 12, kernel_size = 2, stride = 2, padding = 0), 6)
  expect_equal(res_calculate("pooling", in_res = 12, kernel_size = 2), 6)
  expect_equal(res_calculate("layer", in_res = 11, kernel_size = 3, stride = 1, padding = 0), 9)
})

test_that("family_selector", {
  fam <- suppressMessages(family_selector(test_sp, "ind_ha"))
  expect_s3_class(fam, "tbl_df")
  expect_true(all(c("family_name", "family_call") %in% names(fam)))
  expect_gt(nrow(fam), 0)
})

test_that("croppin_hood", {
  crop <- croppin_hood(test_sp[1, ], "x", "y", test_envar, size = 2)
  expect_s4_class(crop, "SpatRaster")
  expect_equal(dim(crop), c(5, 5, terra::nlyr(test_envar)))
  expect_error(
    croppin_hood(test_sp[1, ], "x", "y", test_envar, size = 2, raster_padding = TRUE),
    "Padding method"
  )
})

test_that("cnn_make_samples", {
  s <- cnn_make_samples(test_sp[1:3, ], "x", "y", "ind_ha", test_envar, size = 2)
  expect_named(s, c("predictors", "response"))
  expect_length(s$predictors, 3)
  expect_length(s$response, 3)
  expect_equal(dim(s$predictors[[1]]), c(5, 5, terra::nlyr(test_envar)))
})

test_that("get_partition_samples", {
  s <- get_partition_samples(test_sp, "x", "y", "ind_ha",
    folds = 1:2, partition = ".part1", rasters = test_envar, crop_size = 1
  )
  expect_named(s, c("1", "2"))
  expect_equal(length(s[["1"]]$response), sum(test_sp$.part1 == 1))
})

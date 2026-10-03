small_env <- terra::crop(test_envar, terra::ext(-60, -55, -28, -24))

test_that("adm_predict handles all algorithms individually and as a list", {
  ms <- get_test_models()
  for (nm in names(ms)) {
    p <- suppressWarnings(adm_predict(ms[[nm]], small_env, training_data = test_sp))
    expect_equal(names(p), nm)
    expect_s4_class(p[[1]], "SpatRaster")
  }
  p <- suppressWarnings(adm_predict(ms[c("raf", "glm", "svm")], small_env, training_data = test_sp))
  expect_equal(sort(names(p)), c("glm", "raf", "svm"))
})

test_that("adm_predict options: predict_area, nchunk, inverse transform, quantiles", {
  ms <- get_test_models()
  area <- terra::as.polygons(terra::ext(-59, -56, -27, -25), crs = "EPSG:4326")
  p_area <- adm_predict(ms$raf, small_env, training_data = test_sp, predict_area = area)
  p_all <- adm_predict(ms$raf, small_env, training_data = test_sp)
  expect_lt(
    sum(!is.na(terra::values(p_area[[1]]))),
    sum(!is.na(terra::values(p_all[[1]])))
  )

  p_chunk <- adm_predict(ms$raf, small_env, training_data = test_sp, nchunk = 3)
  expect_equal(terra::values(p_chunk[[1]]), terra::values(p_all[[1]]))

  p_inv <- suppressWarnings(adm_predict(ms$raf, small_env,
    training_data = test_sp,
    invert_transform = c(method = "zscore", a = 10, b = 2)
  ))
  expect_s4_class(p_inv[[1]], "SpatRaster")

  p_neg <- adm_predict(ms$glm, small_env, training_data = test_sp, transform_negative = TRUE)
  expect_gte(terra::global(p_neg[[1]], "min", na.rm = TRUE)[[1]], 0)

  skip_if_not_installed("quantregForest")
  p_q <- adm_predict(ms$qrf, small_env, training_data = test_sp, pred_quantile = 0.75)
  expect_equal(names(p_q), "qrf")
})

test_that("adm_predict rejects invalid models", {
  expect_error(adm_predict(list(a = 1), small_env, training_data = test_sp))
})

test_that("adm_uncertainty works for all supported algorithms", {
  ms <- get_test_models()
  for (nm in setdiff(names(ms), "qrf")) {
    u <- suppressMessages(suppressWarnings(
      adm_uncertainty(ms[[nm]], test_sp, "ind_ha", small_env, iteration = 2)
    ))
    expect_s4_class(u, "SpatRaster")
    expect_equal(terra::nlyr(u), 1)
  }
})

test_that("adm_uncertainty rejects unsupported models", {
  skip_if_not_installed("quantregForest")
  ms <- get_test_models()
  expect_error(adm_uncertainty(ms$qrf, test_sp, "ind_ha", small_env, iteration = 2))
})

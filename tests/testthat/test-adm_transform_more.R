dd <- data.frame(v = c(0, 1, 5, 10, 20))

test_that("adm_transform inverts data.frame transformations", {
  z <- adm_transform(dd, "v", "zscore")
  inv <- adm_transform(z, "v_zscore", "zscore",
    inverse = TRUE, t_terms = c(mean(dd$v), sd(dd$v))
  )
  expect_equal(inv$v_zscore_inverted, dd$v)

  r <- adm_transform(dd, "v", "01")
  inv <- adm_transform(r, "v_01", "01", inverse = TRUE, t_terms = c(min(dd$v), max(dd$v)))
  expect_equal(inv$v_01_inverted, dd$v)

  l <- adm_transform(dd, "v", "log1")
  expect_warning(inv <- adm_transform(l, "v_log1", "log1", inverse = TRUE))
  expect_equal(inv$v_log1_inverted, dd$v)

  l <- adm_transform(dd[-1, , drop = FALSE], "v", "log")
  expect_warning(inv <- adm_transform(l, "v_log", "log", inverse = TRUE))
  expect_equal(inv$v_log_inverted, dd$v[-1])
})

test_that("adm_transform data.frame errors", {
  expect_error(adm_transform(dd, "v", "zscore", inverse = TRUE), "terms")
  expect_error(adm_transform(dd, "v", "round", inverse = TRUE), "cannot be inverted")
  expect_error(adm_transform(dd, "v", "nomethod"), "Undefined")
  expect_error(adm_transform(dd, "v", "nomethod", inverse = TRUE, t_terms = c(1, 1)), "Undefined")
})

test_that("adm_transform works with SpatRaster", {
  r <- terra::rast(nrows = 4, ncols = 4, vals = 1:16)
  names(r) <- "a"
  for (m in c("01", "zscore", "log", "log1", "round")) {
    out <- adm_transform(r, "a", m)
    expect_true(paste0("a_", m) %in% names(out))
  }
  expect_error(adm_transform(r, "a", "nomethod"), "Undefined")

  inv <- adm_transform(adm_transform(r, "a", "zscore"), "a_zscore", "zscore",
    inverse = TRUE, t_terms = c(mean(1:16), sd(1:16))
  )
  expect_equal(terra::values(inv[["a_zscore_inverted"]])[, 1], 1:16, ignore_attr = TRUE)

  inv <- adm_transform(adm_transform(r, "a", "01"), "a_01", "01",
    inverse = TRUE, t_terms = c(1, 16)
  )
  expect_equal(terra::values(inv[["a_01_inverted"]])[, 1], 1:16, ignore_attr = TRUE)

  expect_warning(inv <- adm_transform(adm_transform(r, "a", "log"), "a_log", "log", inverse = TRUE))
  expect_warning(inv <- adm_transform(adm_transform(r, "a", "log1"), "a_log1", "log1", inverse = TRUE))
  expect_error(adm_transform(r, "a", "zscore", inverse = TRUE), "terms")
  expect_error(adm_transform(r, "a", "round", inverse = TRUE), "cannot be inverted")
})

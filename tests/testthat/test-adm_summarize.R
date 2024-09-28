test_that("adm_summarize", {
  require(dplyr)
  require(gamlss)
  envar <- system.file("external/envar.tif", package = "adm")
  envar <- terra::rast(envar)
  data("sppabund")

  # Extract data
  some_sp <- sppabund %>%
    dplyr::filter(species == "Species one") %>%
    dplyr::select(species, ind_ha, x, y, .part)

  # Species abundance data, coordinates, and partition
  some_sp <- sppabund %>%
    dplyr::filter(species == "Species one") %>%
    dplyr::select(species, ind_ha, x, y, .part)
  some_sp

  # Extract data
  some_sp <-
    adm_extract(
      data = some_sp,
      x = "x",
      y = "y",
      env_layer = envar
    )

  # Fit RA
  m_raf <- fit_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("elevation", "sand", "bio3", "bio12"),
    partition = ".part",
    distribution = NO()
  )

  # Fit SVM
  m_svm <- fit_abund_svm(
    data = some_sp,
    response = "ind_ha",
    predictors = c("elevation", "sand", "bio3", "bio12"),
    partition = ".part"
  )

  # XGB
  m_xbg <- fit_abund_xgb(
    data = some_sp,
    response = "ind_ha",
    predictors = c("elevation", "sand", "bio3", "bio12"),
    partition = ".part"
  )

  perf <- adm_summarize(list(m_svm, m_raf, m_xbg))

  # Test the function returns a data frame with a model_ID column when the input contains only one list
  expect_true("model_ID" %in% names(perf))

  # Test the function returns a data frame with sane dimensions
  expect_equal(c(3, 14), dim(perf))

  perf <- adm_summarize(list(m_svm))

  # Test the function returns a data frame with sane dimensions
  expect_equal(c(1, 14), dim(perf))
})

test_that("multiplication works", {
  expect_error(adm_summarize(1))
})

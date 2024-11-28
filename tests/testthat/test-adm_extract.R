test_that("adm_extract performs as expected", {
  require(terra)

  # Load datasets
  data("sppabund")
  envar <- system.file("external/envar.tif", package = "adm")
  envar <- terra::rast(envar)

  # Extract data
  some_sp <- sppabund %>%
    filter(species == "Species one") %>%
    dplyr::select(species, ind_ha, x, y)

  # Run the function
  result <-
    adm_extract(
      data = some_sp,
      x = "x",
      y = "y",
      env_layer = envar,
      variables = NULL,
      filter_na = FALSE
    )


  # Compare the result with the expected output
  expect_equal(dim(result), c(938, 12))

  # Add other tests needed, for example testing that the output is a data.frame
  expect_equal(class(result)[1], "tbl_df")

  # Also, you may wish to test that the result has the correct column names
  expect_equal(
    colnames(result),
    c(
      "species",
      "ind_ha",
      "x",
      "y",
      "bio1",
      "bio12",
      "bio15",
      "bio3",
      "cfvo",
      "elevation",
      "sand",
      "eco"
    )
  )
})

test_that("adm_extract remove NA with lines", {
  require(terra)

  # Load datasets
  data("sppabund")
  envar <- system.file("external/envar.tif", package = "adm")
  envar <- terra::rast(envar)

  # Extract data
  some_sp <- sppabund %>%
    filter(species == "Species one") %>%
    dplyr::select(species, ind_ha, x, y)

  some_sp[1, c("x", "y")] <- some_sp[1, c("x", "y")] * -1

  # Run the function
  result <-
    adm_extract(
      data = some_sp,
      x = "x",
      y = "y",
      env_layer = envar,
      variables = NULL,
      filter_na = TRUE
    )

  expect_message(adm_extract(
    data = some_sp,
    x = "x",
    y = "y",
    env_layer = envar,
    variables = NULL,
    filter_na = TRUE
  ))

  expect_equal(dim(result), c(937, 12))
})

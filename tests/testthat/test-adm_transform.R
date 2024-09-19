data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species one") %>%
  dplyr::select(species, ind_ha, x, y)

test_that("Test for transformation type 01", {
  
  result <- adm_transform(some_sp, variable = "ind_ha", method = "01")
  expect_equal(range(result$ind_ha_01), c(0, 1))
  
})

test_that("Test for transformation type ln", {
  
  result <- adm_transform(some_sp, variable = "ind_ha", method = "log") 
  expect_equal(result$ind_ha_log, log(result$ind_ha))
  
})

test_that("Test for transformation type round", {
  
  result <- adm_transform(some_sp, variable = "ind_ha", method = "round") 
  expect_equal(result$ind_ha_round, round(result$ind_ha))
  
})

test_that("Test for transformation type zscore", {
  
  result <- adm_transform(some_sp, variable = "ind_ha", method = "zscore") 
  expect_equal(round(mean(result$ind_ha_zscore)),  16)
  expect_equal(round(sd(result$ind_ha_zscore)),  32)
})


# Clean up
rm(some_sp)



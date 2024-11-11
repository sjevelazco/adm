require(dplyr)

data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species three") %>%
  dplyr::select(species, ind_ha, x, y)

test_that("absence_ratio 0.5", {
  some_sp_2 <- balance_dataset(
    data = some_sp,
    response = "ind_ha",
    absence_ratio = 0.5
  )
  
  expect_equal(
    nrow(some_sp_2),
    672
  )
})

test_that("absence_ratio 0", {
  some_sp_2 <- balance_dataset(
    data = some_sp,
    response = "ind_ha",
    absence_ratio = 0
  )
  
  expect_equal(
    nrow(some_sp_2),
    448
  )
})

test_that("absence_ratio 0", {
  some_sp_2 <- balance_dataset(
    data = some_sp,
    response = "ind_ha",
    absence_ratio = 0
  )
  
  expect_message(
    balance_dataset(
      data = some_sp_2,
      response = "ind_ha",
      absence_ratio = 0
    )
  )
  expect_message(
    balance_dataset(
      data = some_sp_2,
      response = "ind_ha",
      absence_ratio = 0.5
    )
  )
})

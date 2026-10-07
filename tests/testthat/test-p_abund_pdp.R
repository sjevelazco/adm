test_that("data_abund_pdp and p_abund_pdp", {
  m <- fit_abund_raf(test_sp, "ind_ha", test_pred,
    partition = ".part1", ntree = 20, verbose = FALSE
  )
  d <- data_abund_pdp(m, "bio12",
    resolution = 8,
    training_data = test_sp, projection_data = test_envar
  )
  expect_true("pdpdata" %in% names(d))
  expect_true(all(c("bio12", "Abundance", "Type") %in% names(d$pdpdata)))

  expect_error(data_abund_pdp(m, c("bio12", "sand"),
    resolution = 8,
    training_data = test_sp, projection_data = test_envar
  ))

  p <- p_abund_pdp(m, test_pred[1:2],
    resolution = 8, resid = TRUE,
    training_data = test_sp, projection_data = test_envar
  )
  expect_s3_class(p, "ggplot")
  expect_error(p_abund_pdp(m$model, test_pred[1],
    training_data = test_sp,
    projection_data = test_envar
  ))
})

test_that("p_abund_bpdp", {
  m <- fit_abund_raf(test_sp, "ind_ha", test_pred,
    partition = ".part1", ntree = 20, verbose = FALSE
  )
  p <- p_abund_bpdp(m, c("bio12", "sand"),
    resolution = 8,
    training_data = test_sp, projection_data = test_envar
  )
  expect_s3_class(p, "ggplot")
})

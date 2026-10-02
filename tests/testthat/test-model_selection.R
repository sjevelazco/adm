test_that("model_selection returns one optimal combination", {
  set.seed(1)
  hc <- data.frame(
    comb_id = paste0("comb_", 1:4),
    mtry = c(2, 3, 2, 3),
    model = "raf",
    mae_mean = c(21, 22, 20, 21.5),
    mae_sd = 1,
    corr_spear_mean = c(0.25, 0.22, 0.30, 0.24),
    corr_spear_sd = 0.02,
    corr_pear_mean = c(0.15, 0.12, 0.18, 0.14),
    corr_pear_sd = 0.03,
    inter_mean = c(17, 18, 16, 18),
    inter_sd = 2,
    slope_mean = c(0.3, 0.25, 0.35, 0.27),
    slope_sd = 0.1,
    pdisp_mean = c(0.5, 0.52, 0.5, 0.51),
    pdisp_sd = 0.1
  )
  res <- model_selection(hc, c("corr_spear", "mae"))
  expect_named(res, c("optimal_combination", "all_combinations"))
  expect_equal(nrow(res$optimal_combination), 1)
  expect_equal(res$optimal_combination$comb_id, "comb_3")
  expect_equal(nrow(res$all_combinations), 4)
})

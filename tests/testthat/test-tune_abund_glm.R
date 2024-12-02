# require(dplyr)
# data("sppabund")
# some_sp <- sppabund %>%
#   dplyr::filter(species == "Species two") %>%
#   dplyr::select(-.part2, -.part3)
#
# suitable_distributions <- family_selector(data = some_sp, response = "ind_ha")
#
# # Create a grid
# glm_grid <- expand.grid(
#  poly = c(2, 3),
#  inter_order = c(1, 2),
#  distribution = suitable_distributions[["family_call"]]
# )
#
# test_that("tune_abund_glm", {
#   set.seed(123)
#   tuned_ <- tune_abund_glm(
#     data = some_sp,
#     response = "ind_ha",
#     predictors = c("bio12", "elevation", "sand"),
#     predictors_f = c("eco"),
#     partition = ".part",
#     predict_part = TRUE,
#     metrics = c("corr_pear", "mae"),
#     grid = glm_grid,
#     n_cores = 8
#   )
#
#   dim(tuned_$optimal_combination) %>% expect_equal(c(1, 16))
#   expect_true(tuned_$performance$corr_spear_mean < 0.5)
# })
#
# data = some_sp
# response = "ind_ha"
# predictors = c("bio12", "elevation", "sand")
# predictors_f = c("eco")
# fit_formula = NULL
# sigma_formula = ~1
# nu_formula = ~1
# tau_formula = ~1
# partition = ".part"
# predict_part = TRUE
# grid = glm_grid
# metrics = c("corr_pear", "mae")
# n_cores = 1
# verbose = TRUE
#
#
# test_that("test errors", {
#   expect_error(tune_abund_net(
#     data = some_sp,
#     response = "ind_ha",
#     predictors = c("bio12", "elevation", "sand"),
#     predictors_f = c("eco"),
#     partition = ".part",
#     predict_part = TRUE,
#     # metrics = c("corr_pear", "mae"),
#     grid = grid_0,
#     n_cores = 1
#   ))
#
#   expect_error(tune_abund_net(
#     data = some_sp,
#     response = "ind_ha",
#     predictors = c("bio12", "elevation", "sand"),
#     predictors_f = c("eco"),
#     partition = ".part",
#     predict_part = TRUE,
#     metrics = c("corr_pear", "mae"),
#     grid = expand.grid(
#       mtryE = seq(from = 2, to = 3, by = 1),
#       ntreeE = c(100, 500)
#     ),
#     n_cores = 1
#   ))
# })
#
# test_that("incomplete grid", {
#   tuned_ <- tune_abund_net(
#     data = some_sp,
#     response = "ind_ha",
#     predictors = c("bio12", "elevation", "sand"),
#     predictors_f = c("eco"),
#     partition = ".part",
#     predict_part = TRUE,
#     metrics = c("corr_pear", "mae"),
#     grid =  expand.grid(
#       # size = c(4, 8, 12),
#       decay = seq(from = 0, to = 0.4, by = 0.1)
#     ),
#     n_cores = 1
#   )
#
#   expect_true("size" %in% names(tuned_$optimal_combination))
# })

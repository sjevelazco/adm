require(dplyr)
data("sppabund")
some_sp <- sppabund %>%
  dplyr::filter(species == "Species two") %>%
  dplyr::select(-.part2, -.part3)

some_sp <-
  balance_dataset(some_sp, response = "ind_ha", absence_ratio = 0.2)

suitable_distributions <- family_selector(data = some_sp, response = "ind_ha")

# Create a grid
grid_0 <- expand.grid(
  inter = 1:4,
  distribution = c("LO", "NO"), 
  stringsAsFactors = FALSE
)

test_that("tune_abund_gam", {
  set.seed(123)
  tuned_ <- tune_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    fit_formula = formula("ind_ha ~ bio12 + elevation + sand + eco"),
    sigma_formula = formula("ind_ha ~ bio12 + elevation"),
    nu_formula = formula("ind_ha ~ bio12 + elevation"),
    tau_formula = ~1,
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = grid_0,
    n_cores = 3,
    verbose = TRUE
  )
  
  dim(tuned_$optimal_combination) %>% expect_equal(c(1, 16))
  expect_true(tuned_$performance$corr_spear_mean < 0.5)
})


test_that("test errors", {
  expect_error(tune_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    # metrics = c("corr_pear", "mae"),
    grid = grid_0,
    n_cores = 1
  ))
  
  expect_error(tune_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      mtryE = seq(from = 2, to = 3, by = 1),
      ntreeE = c(100, 500)
    ),
    n_cores = 1
  ))
})

test_that("incomplete grid", {
  tuned_ <- tune_abund_gam(
    data = some_sp,
    response = "ind_ha",
    predictors = c("bio12", "elevation", "sand"),
    fit_formula = formula("ind_ha ~ bio12 + elevation + sand + eco"),
    sigma_formula = formula("ind_ha ~ bio12 + elevation"),
    nu_formula = formula("ind_ha ~ bio12 + elevation"),
    tau_formula = ~1,
    predictors_f = c("eco"),
    partition = ".part",
    predict_part = TRUE,
    metrics = c("corr_pear", "mae"),
    grid = expand.grid(
      # inter = 1:4,
      distribution = c("LO", "NO"), 
      stringsAsFactors = FALSE
    ),
    n_cores = 3,
    verbose = TRUE
  )
  
  expect_true("inter" %in% names(tuned_$optimal_combination))
})

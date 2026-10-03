# Extra paths of fit_abund_*: no partition, hold-out sets, custom formulas

fit_args <- list(
  raf = list(ntree = 20),
  glm = list(distribution = "NO"),
  gam = list(distribution = "NO"),
  net = list(size = 3),
  svm = list(),
  gbm = list(distribution = "gaussian", n.trees = 20),
  xgb = list(nrounds = 10)
)

fit_one <- function(algo, ...) {
  args <- c(
    list(
      data = test_sp, response = "ind_ha", predictors = test_pred,
      verbose = FALSE
    ),
    fit_args[[algo]], list(...)
  )
  suppressMessages(suppressWarnings(do.call(paste0("fit_abund_", algo), args)))
}

test_that("fit_abund_* without partition return model, predictors and metadata", {
  for (algo in names(fit_args)) {
    m <- fit_one(algo, partition = NULL)
    expect_true(all(c("model", "predictors", "metadata") %in% names(m)), info = algo)
    expect_equal(m$predictors$model, algo)
    expect_equal(m$metadata$source_function, paste0("fit_abund_", algo))
  }
})

test_that("fit_abund_* work with a hold-out set and categorical predictors", {
  for (algo in setdiff(names(fit_args), "xgb")) {
    m <- fit_one(algo,
      predictors_f = "eco", partition = ".part1",
      predict_part = TRUE, hold_out_set = test_sp[1:30, ]
    )
    expect_true("performance_part_hold_out" %in% names(m), info = algo)
    expect_false(is.null(m$predicted_part), info = algo)
  }
  m <- fit_one("xgb", partition = ".part1", predict_part = TRUE, hold_out_set = test_sp[1:30, ])
  expect_true("performance_part_hold_out" %in% names(m))
})

test_that("fit_abund_* use a custom formula", {
  f <- ind_ha ~ bio12 + sand
  for (algo in c("raf", "svm", "net", "gbm")) {
    m <- fit_one(algo, fit_formula = f, partition = ".part1")
    expect_equal(nrow(m$performance), 1, info = algo)
  }
})

test_that("fit_abund_glm accepts polynomials and interactions", {
  m <- fit_one("glm", partition = ".part1", poly = 2, inter_order = 2)
  expect_equal(nrow(m$performance), 1)
  expect_s3_class(m$model, "gamlss")
})

test_that("fit_abund_svm sigma and kernel arguments", {
  m <- fit_one("svm", partition = ".part1", sigma = 0.1, C = 2)
  expect_equal(nrow(m$performance), 1)
})

test_that("fit_abund_qrf with hold-out-free paths and categorical predictors", {
  skip_if_not_installed("quantregForest")
  m <- suppressMessages(suppressWarnings(fit_abund_qrf(test_sp, "ind_ha", test_pred, "eco",
    partition = ".part1", train_quantiles = c(0.25, 0.75), eval_quantile = 0.5,
    ntree = 20, verbose = FALSE
  )))
  expect_equal(nrow(m$performance), 1)
})

test_that("tune_abund_qrf", {
  skip_if_not_installed("quantregForest")
  tuned <- suppressMessages(suppressWarnings(tune_abund_qrf(
    data = test_sp, response = "ind_ha", predictors = test_pred,
    partition = ".part1", metrics = c("corr_spear", "mae"),
    grid = expand.grid(mtry = 2:3, ntree = 20, nodesize = 5),
    n_cores = 1, verbose = FALSE
  )))
  expect_true(all(c("optimal_combination", "all_combinations") %in% names(tuned)))
  expect_equal(nrow(tuned$optimal_combination), 1)
  expect_equal(nrow(tuned$all_combinations), 2)
  expect_error(tune_abund_qrf(
    data = test_sp, response = "ind_ha", predictors = test_pred,
    partition = ".part1",
    grid = expand.grid(mtryE = 1, ntreeE = 10), n_cores = 1, verbose = FALSE
  ))
})

# Quick tests for fit_abund_* functions (tiny hyperparameters)

check_fit <- function(m, algo) {
  expect_type(m, "list")
  expect_true(all(c("model", "predictors", "performance", "performance_part", "metadata") %in% names(m)))
  expect_equal(m$performance$model, algo)
  expect_true(all(c("mae_mean", "corr_spear_mean", "pdisp_mean") %in% names(m$performance)))
  expect_equal(nrow(m$performance), 1)
  expect_true(all(is.finite(m$performance$mae_mean)))
}

test_that("fit_abund_raf with and without partition", {
  m <- fit_abund_raf(test_sp, "ind_ha", test_pred, "eco",
    partition = ".part1", predict_part = TRUE, ntree = 20, verbose = FALSE
  )
  check_fit(m, "raf")
  expect_s3_class(m$model, "randomForest")
  expect_false(is.null(m$predicted_part))

  m0 <- fit_abund_raf(test_sp, "ind_ha", test_pred,
    partition = NULL, ntree = 20, verbose = FALSE
  )
  expect_named(m0, c("model", "predictors", "metadata"))
})

test_that("fit_abund_raf with hold-out set", {
  m <- fit_abund_raf(test_sp, "ind_ha", test_pred,
    partition = ".part1", hold_out_set = test_sp[1:20, ],
    ntree = 20, verbose = FALSE
  )
  expect_true("performance_part_hold_out" %in% names(m))
})

test_that("fit_abund_glm and fit_abund_gam with hold-out set", {
  for (fn in list(fit_abund_glm, fit_abund_gam)) {
    m <- suppressWarnings(fn(test_sp, "ind_ha", test_pred,
      partition = ".part1", hold_out_set = test_sp[1:30, ],
      distribution = "NO", verbose = FALSE
    ))
    expect_true("performance_part_hold_out" %in% names(m))
    expect_true(all(is.finite(m$performance_part_hold_out$mae)))
  }
})

test_that("fit_abund_glm", {
  m <- suppressWarnings(fit_abund_glm(test_sp, "ind_ha", test_pred,
    partition = ".part1", distribution = "NO", verbose = FALSE
  ))
  check_fit(m, "glm")
  expect_error(
    fit_abund_glm(test_sp, "ind_ha", test_pred, partition = ".part1", verbose = FALSE),
    "distribution"
  )
})

test_that("fit_abund_gam", {
  m <- suppressWarnings(fit_abund_gam(test_sp, "ind_ha", test_pred,
    partition = ".part1", distribution = "NO", verbose = FALSE
  ))
  check_fit(m, "gam")
  expect_error(
    fit_abund_gam(test_sp, "ind_ha", test_pred, partition = ".part1", verbose = FALSE)
  )
})

test_that("fit_abund_net", {
  m <- suppressWarnings(fit_abund_net(test_sp, "ind_ha", test_pred,
    partition = ".part1", size = 3, verbose = FALSE
  ))
  check_fit(m, "net")
})

test_that("fit_abund_svm", {
  m <- suppressWarnings(fit_abund_svm(test_sp, "ind_ha", test_pred,
    partition = ".part1", verbose = FALSE
  ))
  check_fit(m, "svm")
})

test_that("fit_abund_gbm", {
  m <- suppressMessages(suppressWarnings(fit_abund_gbm(test_sp, "ind_ha", test_pred,
    partition = ".part1", distribution = "gaussian", n.trees = 20, verbose = FALSE
  )))
  check_fit(m, "gbm")
})

test_that("fit_abund_xgb", {
  m <- suppressWarnings(fit_abund_xgb(test_sp, "ind_ha", test_pred,
    partition = ".part1", nrounds = 20, verbose = FALSE
  ))
  check_fit(m, "xgb")
  expect_warning(
    fit_abund_xgb(test_sp, "ind_ha", test_pred, predictors_f = "eco",
      partition = ".part1", nrounds = 5, verbose = FALSE),
    "Categorical"
  )
})

test_that("fit_abund_qrf with and without partition", {
  skip_if_not_installed("quantregForest")
  m <- fit_abund_qrf(test_sp, "ind_ha", test_pred,
    partition = ".part1", predict_part = TRUE, ntree = 20, verbose = FALSE
  )
  check_fit(m, "qrf")
  expect_false(is.null(m$predicted_part))

  m0 <- fit_abund_qrf(test_sp, "ind_ha", test_pred,
    partition = NULL, ntree = 20, verbose = FALSE
  )
  expect_named(m0, c("model", "predictors", "metadata"))
})

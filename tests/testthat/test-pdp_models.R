# PDP/BPDP data and plots across model classes (tiny resolution for speed)

test_that("data_abund_pdp works for all algorithms", {
  ms <- get_test_models()
  ms$qrf <- NULL # qrf is not supported by the pdp functions
  for (nm in names(ms)) {
    d <- data_abund_pdp(ms[[nm]], "bio12",
      resolution = 5, resid = TRUE,
      training_data = test_sp, projection_data = test_envar
    )
    expect_true(all(c("pdpdata", "resid") %in% names(d)), info = nm)
    expect_true(nrow(d$pdpdata) > 0, info = nm)
    expect_false(is.null(d$resid), info = nm)
  }
})

test_that("data_abund_pdp without projection data and with categorical predictors", {
  ms <- get_test_models()
  d <- data_abund_pdp(ms$raf, "bio12", resolution = 5, training_data = test_sp)
  expect_true(nrow(d$pdpdata) > 0)
  expect_null(d$resid)

  d <- data_abund_pdp(ms$raf, "eco", training_data = test_sp)
  expect_true(nrow(d$pdpdata) > 1)
  d <- data_abund_pdp(ms$glm, "eco", training_data = test_sp)
  expect_true(nrow(d$pdpdata) > 1)
})

test_that("data_abund_pdp inverts transformations and checks inputs", {
  ms <- get_test_models()
  expect_warning(
    d <- data_abund_pdp(ms$raf, "bio12",
      resolution = 5, training_data = test_sp,
      invert_transform = c(method = "log", a = NA, b = NA)
    ),
    "log"
  )
  expect_s3_class(d$pdpdata, "data.frame")

  expect_error(data_abund_pdp(ms$glm, "bio12", resolution = 5), "training_data")
  expect_error(data_abund_pdp(ms$raf, "bio12", resolution = 5), "training_data")
  expect_error(data_abund_pdp(ms$raf$model, "bio12", training_data = test_sp))
})

test_that("data_abund_bpdp works for all algorithms", {
  ms <- get_test_models()
  ms$qrf <- NULL
  for (nm in names(ms)) {
    d <- data_abund_bpdp(ms[[nm]], c("bio12", "sand"),
      resolution = 5,
      training_data = test_sp, projection_data = test_envar,
      training_boundaries = "convexh"
    )
    expect_type(d, "list")
    expect_true(length(d) > 0, info = nm)
  }
})

test_that("data_abund_bpdp boundary types and input checks", {
  ms <- get_test_models()
  for (b in c("rectangle", "convexh", "none")) {
    d <- try(data_abund_bpdp(ms$raf, c("bio12", "sand"),
      resolution = 5,
      training_data = test_sp, projection_data = test_envar,
      training_boundaries = b
    ), silent = TRUE)
    expect_false(inherits(d, "try-error") && b != "none", info = b)
  }
  expect_error(data_abund_bpdp(ms$raf, "bio12", resolution = 5, training_data = test_sp))
  expect_error(data_abund_bpdp(ms$glm, c("bio12", "sand"), resolution = 5))
})

test_that("p_abund_pdp and p_abund_bpdp draw plots for several algorithms", {
  # models without categorical predictors (not available as raster layers)
  ms <- list(
    glm = suppressWarnings(fit_abund_glm(test_sp, "ind_ha", test_pred,
      partition = ".part1", distribution = "NO", verbose = FALSE
    )),
    net = suppressWarnings(fit_abund_net(test_sp, "ind_ha", test_pred,
      partition = ".part1", size = 3, verbose = FALSE
    )),
    xgb = get_test_models()$xgb
  )
  for (nm in c("glm", "net", "xgb")) {
    p <- p_abund_pdp(ms[[nm]], c("bio12", "sand"),
      resolution = 5, resid = (nm == "xgb"), rug = TRUE,
      training_data = test_sp, projection_data = test_envar,
      response_name = "Abundance"
    )
    expect_s3_class(p, "ggplot")
    p <- p_abund_bpdp(ms[[nm]], c("bio12", "sand"),
      resolution = 5,
      training_data = test_sp, projection_data = test_envar,
      response_name = "Abundance"
    )
    expect_s3_class(p, "ggplot")
  }
  ms <- get_test_models()
  # all predictors / all pairs when predictors = NULL
  p <- p_abund_pdp(ms$raf, NULL, resolution = 5, training_data = test_sp)
  expect_s3_class(p, "ggplot")
  p <- p_abund_bpdp(ms$raf, NULL, resolution = 5, training_data = test_sp)
  expect_s3_class(p, "ggplot")
  # categorical predictor
  p <- p_abund_pdp(ms$raf, "eco", training_data = test_sp)
  expect_s3_class(p, "ggplot")
})

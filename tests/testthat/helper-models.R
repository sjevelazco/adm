# Tiny fitted models shared by several test files (built once, lazily)
.model_cache <- new.env()

get_test_models <- function() {
  if (!is.null(.model_cache$models)) {
    return(.model_cache$models)
  }
  p <- ".part1"
  quiet <- function(x) suppressMessages(suppressWarnings(x))
  m <- list(
    raf = quiet(fit_abund_raf(test_sp, "ind_ha", test_pred, "eco",
      partition = p, ntree = 20, verbose = FALSE
    )),
    glm = quiet(fit_abund_glm(test_sp, "ind_ha", test_pred, "eco",
      partition = p, distribution = "NO", poly = 2, verbose = FALSE
    )),
    gam = quiet(fit_abund_gam(test_sp, "ind_ha", test_pred, "eco",
      partition = p, distribution = "NO", verbose = FALSE
    )),
    net = quiet(fit_abund_net(test_sp, "ind_ha", test_pred, "eco",
      partition = p, size = 3, verbose = FALSE
    )),
    svm = quiet(fit_abund_svm(test_sp, "ind_ha", test_pred, "eco",
      partition = p, verbose = FALSE
    )),
    gbm = quiet(fit_abund_gbm(test_sp, "ind_ha", test_pred, "eco",
      partition = p, distribution = "gaussian", n.trees = 20, verbose = FALSE
    )),
    xgb = quiet(fit_abund_xgb(test_sp, "ind_ha", test_pred,
      partition = p, nrounds = 20, verbose = FALSE
    ))
  )
  if (requireNamespace("quantregForest", quietly = TRUE)) {
    m$qrf <- quiet(fit_abund_qrf(test_sp, "ind_ha", test_pred, "eco",
      partition = p, ntree = 20, verbose = FALSE
    ))
  }
  .model_cache$models <- m
  m
}

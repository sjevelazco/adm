#' Merge model performance tables
#'
#' @param models list. A list of a single or several models fitted with some of fit_ or tune_ functions. Usage models = list(mod1, mod2, mod3)
#'
#' @importFrom dplyr bind_rows relocate tibble select %>%
#'
#' @return A tibble object with combined model performance for all input models. Models fit with tune will include model performance for the best hyperparameters.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#' 
#' data("sppabund")
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)
#' 
#' # Species abundance data, coordinates, and partition
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(species, ind_ha, x, y, .part)
#' some_sp
#' 
#' # Extract data
#' some_sp <-
#'   adm_extract(
#'     data = some_sp,
#'     x = "x",
#'     y = "y",
#'     env_layer = envar
#'   )
#' 
#' # Fit RAF
#' m_raf <- fit_abund_raf(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   partition = ".part",
#' )
#' 
#' # Fit SVM
#' m_svm <- fit_abund_svm(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   partition = ".part"
#' )
#' 
#' # XGB
#' m_xbg <- fit_abund_xgb(
#'   data = some_sp,
#'   response = "ind_ha",
#'   predictors = c("elevation", "sand", "bio3", "bio12"),
#'   partition = ".part"
#' )
#' 
#' 
#' perf <- adm_summarize(list(m_svm, m_raf, m_xbg))
#' 
#' perf
#' }
adm_summarize <- function(models) {
  . <- model_ID <- model <- pdisp_sd <- NULL
  if (data.class(models) != "list") {
    stop("models must be a list object")
  }
  if (length(models) > 1) {
    perf <- lapply(models, function(x) {
      x$performance
    })
    perf <- Map(cbind, perf, model_ID = 1:length(perf))
    perf_tib <- dplyr::bind_rows(perf) %>%
      dplyr::relocate(model_ID, .before = model) %>%
      dplyr::tibble()
  } else {
    perf_tib <- models[[1]]$performance
    perf_tib$model_ID <- 1
  }
  perf_tib <- perf_tib %>% dplyr::relocate(names(dplyr::select(
    perf_tib,
    model_ID:pdisp_sd
  )))
  return(perf_tib)
}


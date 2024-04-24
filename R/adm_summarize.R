#' Merge model performance tables
#'
#' @param models list. A list a single or several models fitted with some of fit_ or tune_ functions. Usage models = list(mod1, mod2, mod3)
#'
#' @importFrom dplyr bind_rows relocate tibble select %>%
#'
#' @return Combined model performance table for all input models. Models fit with tune will include model performance for the best hyperparameters.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' require(flexsdm)
#' data("sppabund")
#' 
#' sp1 <- sppabund %>% 
#'   dplyr::filter(species == "Species one", ind_ha>0) %>% 
#'   dplyr::mutate(ind_ha = ind_ha %>% round())
#' # Families for GAM and GLM
#' dis_f <- system.file("external/families_bank.txt", package = "adm") %>%
#'   utils::read.delim(., header = TRUE, quote = "\t") %>% 
#'   dplyr::as_tibble()
#' dis_f$family_name
#' 
#' # Fit GAM 
#' m_gam <- fit_abund_gam(
#'   data = sp1,
#'   response = "ind_ha",
#'   predictors = c("bio1", "bio12", "bio15", "bio3", "cfvo", "elevation"),
#'   partition = ".part",
#'   family = "PO", 
#'   inter = 1
#' )
#' 
#' m_svm <- fit_abund_svm(
#'   data = sp1,
#'   response = "ind_ha",
#'   predictors = c("bio1", "bio12", "bio15", "bio3", "cfvo", "elevation"),
#'   partition = ".part"
#' )
#' 
#' m_raf <- fit_abund_raf(
#'   data = sp1,
#'   response = "ind_ha",
#'   predictors = c("bio1", "bio12", "bio15", "bio3", "cfvo", "elevation"),
#'   partition = ".part"
#' )
#' 
#' adm_summarize(list(m_gam, m_svm, m_raf))
#' 
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

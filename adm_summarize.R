adm_summarize <- function (models) 
{
  . <- model_ID <- model <- pdispersion_sd <- NULL
  if (data.class(models) != "list") {
    stop("models must be a list object")
  }
  if (length(models) > 1) {
    perf <- lapply(models, function(x) {
      x$performance
    })
    perf <- Map(cbind, perf, model_ID = 1:length(perf))
    perf_tib <- dplyr::bind_rows(perf) %>% dplyr::relocate(model_ID, 
                                                           .before = model) %>% dplyr::tibble()
  } else {
    perf_tib <- models[[1]]$performance
    perf_tib$model_ID <- 1
  }
  perf_tib <- perf_tib %>% dplyr::relocate(names(dplyr::select(perf_tib, 
                                                               model_ID:pdispersion_sd)))
  return(perf_tib)
}

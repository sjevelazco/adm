#' Title TODO
#'
#' @param hyper_combinations 
#' @param metrics 
#'
#' @importFrom dplyr filter select
#'
#' @return
#' @export
#'
#' @examples
model_selection <- function(hyper_combinations, metrics) {
  metrics <- as.vector(metrics)
  performance_var <- names(hyper_combinations)

  mae_columns <- grep("^mae_mean", performance_var)
  hyper_combinations[, mae_columns] <- -1 * hyper_combinations[, mae_columns]

  performance_dict <- list(
    corr_spear = grep("^corr_spear", performance_var),
    corr_pear = grep("^corr_pear", performance_var),
    mae = grep("^mae", performance_var),
    inter = grep("^inter", performance_var),
    slope = grep("^slope", performance_var),
    pdisp = grep("^pdisp", performance_var)
  )
  
  cols_idx <- unlist(performance_dict[metrics])
  perf_means <- performance_var[cols_idx[cols_idx %in% grep("_mean$", performance_var)]]
  if ("pdisp" %in% metrics){
    hyper_combinations$pdisp_dist <- ((1-hyper_combinations$pdisp_mean %>% abs())*-1)
    perf_means[which(perf_means=="pdisp_mean")] <- "pdisp_dist"
  }
  
  not_selected_combinations <- hyper_combinations
  hyper_combinations <- hyper_combinations[!duplicated(hyper_combinations%>%dplyr::select(-comb_id)),]
  while (nrow(hyper_combinations) > 1) {
    for (i in perf_means) {
      hyper_combinations <- hyper_combinations[which(hyper_combinations[i] >= summary(hyper_combinations[[i]])[5]), ]
      if (nrow(hyper_combinations == 1)) {
        break
      }
    }
  }

  selected_comb <- hyper_combinations$comb_id[1]

  # not_selected_combinations <- not_selected_combinations %>%
  #   dplyr::filter(comb_id != selected_comb)

  hyper_combinations[, mae_columns] <- -1 * hyper_combinations[, mae_columns]  
  not_selected_combinations[, mae_columns] <- -1 * not_selected_combinations[, mae_columns]
  
  if ("pdisp" %in% metrics){
    hyper_combinations <- hyper_combinations %>% 
      dplyr::select(-pdisp_dist)
    not_selected_combinations <- not_selected_combinations %>%
      dplyr::select(-pdisp_dist)
  }
  
  return_list <- list(
    "optimal_combination" = hyper_combinations,
    "all_combinations" = not_selected_combinations
  )

  return(return_list)
}

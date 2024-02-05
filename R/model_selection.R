model_selection <- function(hyper_combinations, metrics) {
  metrics <- as.vector(metrics)
  performance_var <- names(hyper_combinations)

  mae_columns <- grep("^mae_mean", performance_var)
  hyper_combinations[, mae_columns] <- -1 * hyper_combinations[, mae_columns]

  performance_dict <- list(
    spearman = grep("^corr_spear", performance_var),
    pearson = grep("^corr_pear", performance_var),
    mae = grep("^mae", performance_var),
    intercept = grep("^inter", performance_var),
    slope = grep("^slope", performance_var),
    p_dispersion = grep("^pdisp", performance_var)
  )

  cols_idx <- unlist(performance_dict[metrics])
  perf_means <- performance_var[cols_idx[cols_idx %in% grep("_mean$", performance_var)]]

  not_selected_combinations <- hyper_combinations
  while (nrow(hyper_combinations) > 1) {
    for (i in perf_means) {
      hyper_combinations <- hyper_combinations[which(hyper_combinations["corr_spear_mean"] >= summary(hyper_combinations[["corr_spear_mean"]])[5]), ]
      if (nrow(hyper_combinations == 1)) {
        break
      }
    }
  }

  selected_comb <- hyper_combinations$comb_id[1]

  not_selected_combinations <- not_selected_combinations %>%
    dplyr::filter(comb_id != selected_comb)

  return_list <- list(
    "optimal_combination" = hyper_combinations,
    "not_optimal_combinations" = not_selected_combinations
  )

  return(return_list)
}

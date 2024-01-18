model_selection <-
  function(hyper_combinations,metrics){
    metrics <- as.vector(metrics)
    performance_var <- names(hyper_combinations)
    
    mae_columns <- grep("^mae_mean",performance_var)
    
    hyper_combinations[,mae_columns] <-  -1*hyper_combinations[,mae_columns]
    
    performance_dict <- list(
      spearman = grep("^corr_spear",performance_var),
      pearson = grep("^corr_pear",performance_var),
      mae = grep("^mae",performance_var),
      intercept = grep("^inter",performance_var),
      slope = grep("^slope",performance_var),
      p_dispersion = grep("^pdispersion",performance_var)
    )
  
    cols_idx <- unlist(performance_dict[metrics])
    means <- performance_var[cols_idx[cols_idx %in% grep("_mean$",performance_var)]]
    sds <- performance_var[cols_idx[cols_idx %in% grep("sd$",performance_var)]]
    
    ranked_combinations <- hyper_combinations %>% 
      dplyr::arrange(
        dplyr::across(
          dplyr::all_of(means), desc),
        dplyr::across(
          dplyr::all_of(sds))
        )

    return(ranked_combinations)
  }

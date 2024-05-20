#' Fit and validate Generalized Linear Models with exploration of hyper-parameters that optimize performance
#'
#' @param data 
#' @param response 
#' @param predictors 
#' @param predictors_f 
#' @param fit_formula 
#' @param partition 
#' @param predict_part 
#' @param grid 
#' @param metrics 
#' @param n_cores 
#' @param verbose 
#'
#' @importFrom doSNOW registerDoSNOW
#' @importFrom dplyr bind_rows left_join select
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats na.omit
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar
#' 
#' @return
#' @export
#'
#' @examples
tune_abund_glm <-
  function(data,
           response,
           predictors,
           predictors_f = NULL,
           fit_formula = NULL,
           partition,
           predict_part = FALSE,
           grid = NULL,
           metrics = NULL,
           n_cores = 1,
           verbose = FALSE) {
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f)) %>%
      t() %>%
      as.vector()
    
    families_bank <-
      system.file("external/families_bank.txt", package = "adm") %>%
      utils::read.delim(., header = TRUE, quote = "\t") # families_bank
    
    # making grid
    if (is.null(grid)){
      message("Grid not provided. Using the default one for GLM.")
      families_hp <- family_selector(data,response)
      grid <- list(
        poly = c(1,2,3),
        inter_order = c(0,1,2),
        family_call = families_hp$family_call
      ) %>%
        expand.grid() %>%
        dplyr::left_join(families_hp, by = "family_call")
    } else if (is.data.frame(grid) & all(names(grid)%in%c("family_call","poly","inter_order")) & all(grid$family_call %in% families_bank$family_call)){
      message("Testing with provided grid.")
      grid <- dplyr::left_join(grid,families_bank,by="family_call") %>% 
        dplyr::select(poly,inter_order,family_call,discrete)
    } else {
      stop("Grid names expected to be 'family_call', 'poly' and 'inter_order'.")
    }
    
    comb_id <- paste("comb_", 1:nrow(grid), sep = "")
    grid <- cbind(comb_id,grid)
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    
    cl <- parallel::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    pb <- utils::txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_glm","adm_eval"), .packages = c("dplyr")) %dopar% {
      data_fam <- data
      if (grid[i,"discrete"]==1){
        data_fam[,response] <- round(data[,response])
      }
      
      model <- tryCatch({
        model <-
          fit_abund_glm(
            data = data_fam,
            response = response,
            predictors = predictors,
            predictors_f = predictors_f,
            fit_formula = fit_formula,
            partition = partition,
            predict_part = predict_part,
            family = grid[i,"family_call"],
            poly = grid[i,"poly"],
            inter_order = grid[i,"inter_order"],
            verbose = FALSE
          )
      }, error = function(err) {
        print("error")
        model <- list(performance = "error")
      })
      
      l <- list(cbind(grid[i,c("comb_id","family_call","poly","inter_order")], model$performance))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)
    
    hyper_combinations <- lapply(hyper_combinations, function(x) dplyr::bind_rows(x)) %>% 
      dplyr::bind_rows()
    
    if ("model$performance" %in% names(hyper_combinations)){
      hyper_combinations <- hyper_combinations %>% 
        dplyr::select(-`model$performance`) 
    }
    
    hyper_combinations <- hyper_combinations %>% 
      stats::na.omit()
    
    row.names(hyper_combinations) <- NULL
    
    ranked_combinations <- model_selection(hyper_combinations, metrics)
    
    # fit final model
    
    choosen_family <- ranked_combinations[[1]][1,"family_call"]
    choosen_poly <- ranked_combinations[[1]][1,"poly"]
    choosen_inter_order <- ranked_combinations[[1]][1,"inter_order"]
    full_data <- data
    if (families_bank[which(families_bank$family_call == choosen_family),"discrete"]==1){
      full_data[,"ind_ha"] <- round(full_data[,"ind_ha"])
    }
    
    message("Fitting the best model...")
    final_model <-
      fit_abund_glm(
        data = full_data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        family = choosen_family,
        poly = choosen_poly,
        inter_order = choosen_inter_order
      )

    message(
      "The best model was a GLM with:", 
      "\n family = ",
      choosen_family,
      "\n poly = ",
      choosen_poly,
      "\n inter_order = ",
      choosen_inter_order
    )
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }

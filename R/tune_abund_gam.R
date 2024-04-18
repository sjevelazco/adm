#' Fit and validate Generalized Additive Models with exploration of hyper-parameters that optimize performance
#'
#' @param data 
#' @param response character. Column name with species abundance data (e.g., 0, 1, 45).
#' @param predictors character. Vector with the column names of quantitative predictor variables
#' (i.e. continuous variables). Usage predictors = c("temperature", "sand", "elevation")
#' @param predictors_f character. Vector with the column names of qualitative predictor variables
#' (i.e. ordinal or nominal variables type). Usage predictors_f = c("landform")
#' @param fit_formula formula. A formula object with response and predictor
#' variables (e.g. formula(abund ~ temp + sand + pH + landform)).
#' Note that the variables used here must be consistent with those used in
#' response, predictors, and predictors_f arguments. Default NULL
#' @param partition  character. Column name with training and validation partition groups.
#' @param predict_part 
#' @param grid 
#' @param metrics 
#' @param n_cores 
#' @param verbose 
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr bind_rows left_join select
#' @importFrom foreach foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats na.omit
#' 
#' @return
#' @export
#'
#' @examples
tune_abund_gam <-
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
    discrete <- i <- NULL
    
    if (is.null(metrics) |
        !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
      stop("Metrics is needed to be defined in 'metric' argument")
    }
    
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f)) %>%
      t() %>%
      as.vector()
    
    data("gamlss_families_table", envir = environment()) # families_bank
    
    # making grid
    if (is.null(grid)){
      message("Families not provided. Picking recommended families to fit data.")
      family_call <- family_selector(data,response)$family_call
      inter <- seq(1,2*length(variables),1)
      grid <- expand.grid(family_call = family_call, inter = inter)
      grid <- dplyr::left_join(grid,families_bank,by="family_call") %>% 
        dplyr::select(family_call,discrete,inter)
    } else if (is.data.frame(grid) & all(names(grid) %in% c("family_call","inter")) & all(grid$family_call %in% families_bank$family_call)) {
      message("Testing with provided families.")
      grid <- dplyr::left_join(grid,families_bank,by="family_call") %>% 
        select(family_call,discrete,inter)
    } else {
      stop("Grid expected to be a vector of gamlss family calls.")
    }
    
    comb_id <- paste("comb_", 1:dplyr::select(grid), sep = "")
    grid <- cbind(comb_id,grid)
    
    # looping the grid
    message("Searching for optimal hyperparameters...")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = nrow(grid), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    hyper_combinations <- foreach::foreach(i = 1:nrow(grid), .options.snow = opts, .export = c("fit_abund_gam","adm_eval"), .packages = c("dplyr","gamlss")) %dopar% {
      data_fam <- data
      if (grid[i,"discrete"]==1){
        data_fam[,response] <- round(data[,response])
      }
      
      model <- tryCatch({
        model <-
          fit_abund_gam(
            data = data,
            response = response,
            predictors = predictors,
            predictors_f = predictors_f,
            fit_formula = fit_formula,
            partition = partition,
            predict_part = predict_part,
            family = grid[i,"family_call"],
            inter = grid[i,"inter"]
          )
      }, error = function(err) {
        print("error")
        model <- list(performance = "error")
      })
      
      l <- list(cbind(grid[i,c("comb_id","family_call","inter")], model[, "performance"]))
      names(l) <- grid[i, "comb_id"]
      l
    }
    parallel::stopCluster(cl)
    
    hyper_combinations <- lapply(hyper_combinations, function(x) bind_rows(x)) %>% 
      dplyr::bind_rows()
    
    hyper_combinations <- hyper_combinations %>% 
      dplyr::select(-`model$performance`) %>% 
      stats::na.omit()
    
    row.names(hyper_combinations) <- NULL
    
    ranked_combinations <- dplyr::bind_rows(hyper_combinations, metrics)
    
    # fit final model
    
    choosen_family <- ranked_combinations[[1]][1,"family_call"]
    full_data <- data
    if (families_bank[which(families_bank$family_call == choosen_family),"discrete"]==1){
      full_data[,"ind_ha"] <- round(full_data[,"ind_ha"])
    }
    
    message("Fitting the best model...")
    final_model <-
      fit_abund_gam(
        data = full_data,
        response = response,
        predictors = predictors,
        predictors_f = predictors_f,
        fit_formula = fit_formula,
        partition = partition,
        predict_part = predict_part,
        family = choosen_family,
        inter = ranked_combinations[[1]][1,"inter"]
      )
    
    message(
      "The best model was a GAM with:", 
      "\n family = ",
      choosen_family,
      "\n inter = ",
      ranked_combinations[[1]][1,"inter"]
    )
    
    final_list <- c(final_model, ranked_combinations)
    
    return(final_list)
  }

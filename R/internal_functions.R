## %######################################################%##
#                                                          #
####                Internal functions                  ####
#                                                          #
## %######################################################%##

#' adapt_df
#'
#' @noRd
adapt_df <- function(data, predictors, predictors_f, response, partition, xy = NULL) {
  data <- data.frame(data)
  if (is.vector(xy)) {
    xy_cols <- data %>%
      dplyr::select(dplyr::all_of(xy))
    xy_cols <- data.frame(xy_cols)
  }
  
  if (is.null(predictors_f)) {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::starts_with(partition))
    data <- data.frame(data)
  } else {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), dplyr::starts_with(partition))
    data <- data.frame(data)
    for (i in predictors_f) {
      data[, i] <- as.factor(data[, i])
    }
  }
  
  if (is.vector(xy)) {
    data <- dplyr::bind_cols(data, xy_cols)
  }
  
  return(data)
}

#' cnn_get_crop_size
#'
#' @noRd
cnn_get_crop_size <- function(sample_size) {
  if (!is.vector(sample_size)) {
    stop("Please, provide a vector containing the sample size c(width,height)")
  } else if (!(sample_size[[1]] == sample_size[[2]])) {
    stop("adm currently only accepts square samples.")
  } else {
    crop_size <- floor(sample_size[[1]] / 2)
  }
  
  return(crop_size)
}

#' get_variables
#'
#' @noRd
get_variables <- function(predictors, predictors_f){
  if (!is.null(predictors_f)) {
    variables <- dplyr::bind_rows(c(c = predictors, f = predictors_f))
  } else {
    variables <- dplyr::bind_rows(c(c = predictors))
  }
  
  return(variables)
}

#' infer_formula
#'
#' @noRd
infer_formula <- function(fit_formula, response, predictors, predictors_f, verbose){
  if (is.null(fit_formula)) {
    formula1 <- stats::formula(paste(response, "~", paste(c(
      predictors,
      predictors_f
    ), collapse = " + ")))
  } else {
    formula1 <- fit_formula
  }
  
  if (verbose) {
    message(
      "Formula used for model fitting:\n",
      Reduce(paste, deparse(formula1)) %>% gsub(paste("  ", "   ", collapse = "|"), " ", .),
      "\n"
    )
  }
  
  return(formula1)
}

#' wrap_final_list
#'
#' @noRd
wrap_final_list <- function(algo, full_model, variables, response, eval_partial_list, predict_part, part_pred_list, metadata){
  # bind predicted evaluation
  eval_partial <- eval_partial_list %>%
    dplyr::bind_rows(.id = "replica") %>%
    dplyr::as_tibble()
  
  # bind predicted partition
  if (predict_part) {
    part_pred <- part_pred_list %>%
      dplyr::bind_rows(.id = "replica")
  } else {
    part_pred <- NULL
  }
  
  # Summarize performance
  eval_final <- eval_partial %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(dplyr::across(mae:pdisp, list(
      mean = mean,
      sd = stats::sd
    )), .groups = "drop")
  
  variables <- dplyr::bind_cols(
    data.frame(
      model = algo,
      response = response
    ),
    variables
  ) %>% as_tibble()
  
  # Final object
  data_list <- list(
    model = full_model,
    predictors = variables,
    performance = eval_final,
    performance_part = eval_partial,
    predicted_part = part_pred
  )
  
  # Standardize output list
  for (i in 2:length(data_list)) {
    if (!class(data_list[[i]])[1] == "tbl_df") {
      data_list[[i]] <- dplyr::as_tibble(data_list[[i]])
    }
  }
  
  data_list <- append(data_list, list(metadata = metadata))
  
  return(data_list)
}

#' get_metadata
#'
#' @noRd
get_metadata <- function(algo, ...){
  metadata <- list(
    plataform = version$platform,
    R_version = version$version.string,
    run_date = Sys.time()
  )
  
  metadata <- append(metadata, list(...)[[1]])
  
  return(metadata)
}

#' check_metrics
#'
#' @noRd
check_metrics <- function(metrics){
  if (is.null(metrics) |
      !all(metrics %in% c("corr_spear", "corr_pear", "mae", "inter", "slope", "pdisp"))) {
    stop("Metrics is needed to be defined in 'metric' argument")
  }
}

#' build_search_grid
#'
#' @noRd
build_search_grid <- function(grid, grid_dict){
  nms_grid <- names(grid)
  nms_hypers <- names(grid_dict)
  
  if (!all(nms_grid %in% nms_hypers)) {
    stop(
      paste(paste(nms_grid[!nms_grid %in% nms_hypers], collapse = ", "), " is not hyperparameters\n"),
      "Grid expected to be any combination between ", paste(nms_hypers, collapse = ", ")
    )
  }
  
  if (is.null(grid)) {
    message("Grid not provided. Using the default.")
    grid <- expand.grid(grid_dict, stringsAsFactors = FALSE)
  } else if (all(nms_hypers %in% nms_grid)) {
    message("Using provided grid.")
  } else if (any(!nms_hypers %in% nms_grid)) {
    message(
      "Adding default hyperparameter for: ",
      paste(names(grid_dict)[!names(grid_dict) %in% nms_grid], collapse = ", ")
    )
    
    user_hyper <- names(grid)[which(names(grid) %in% names(grid_dict))]
    default_hyper <- names(grid_dict)[which(!names(grid_dict) %in% user_hyper)]
    
    user_list <- grid_dict[default_hyper]
    for (i in user_hyper) {
      l <- grid[[i]] %>%
        unique() %>%
        list()
      names(l) <- i
      user_list <- append(user_list, l)
    }
    
    grid <- expand.grid(user_list, stringsAsFactors = FALSE)
  }
  
  return(grid)
}

#' observer_init
#'
#' @noRd
observer_init <- function(){
  list(
    early_stop = c()
  )
}

#' observer_register
#'
#' @noRd
observer_register <- function(observer, what, how){
  if(what == "early_stop"){
    observer[[what]] <- c(observer[[what]], how)
  }
  
  return(observer)
}

#' early_stop_interpreter
#'
#' @noRd
early_stop_interpreter <- function(early_stopping, observer, nrounds){
  early_value <- switch (early_stopping$fm_strategy[[1]],
                         "mean" = {
                           mean(observer$early_stop)
                         },
                         "median" = {
                           median(observer$early_stop)
                         },
                         "min" = {
                           min(observer$early_stop)
                         },
                         "max" = {
                           max(observer$early_stop)
                         },
                         "none" = {
                           nrounds
                         }
  )
  
  return(round(early_value))
}

#' check_models_validity
#'
#' @param models 
#'
#' @noRd
check_models_validity <- function(models){
  # if list check the requirements
  # 1. it is properly organized in sublists
  # 2. each one containing one algorithm
  # 3. and a predictors table
  if(is.list(models)){
    all_list <- all(lapply(models, function(x){class(x) == "list"}) %>% unlist())
    
    if(all_list){
      has_names <- lapply(models, function(x){
        all(c("model","predictors") %in% names(x))
      }) %>% unlist() %>% all()
    } else {
      has_names <- all(c("model","predictors") %in% names(models))
    }
    
    if(all_list & has_names){
      return(c(TRUE, "list_of_models"))
    } else if (has_names) {
      return(c(TRUE, "individual_models"))
    }
  } 
  
  return(FALSE)
}

#' get_predictor_names
#'
#' @param m_detect 
#'
#' @noRd
get_predictor_names <- function(m_detect, i) {
  m_detect[[i]] %>% select(-all_of(c("model","response"))) %>% as.vector() %>% unlist
}


#' filter_safe_levels
#'
#' @param m_detect 
#' @param pred_df 
#' @param training_data 
#'
#' @noRd
filter_safe_levels <- function(m_detect, pred_df, training_data) {
  f_cols <- m_detect[1, grep("f", names(m_detect))] |> unlist() |> as.character()
  
  if (length(f_cols) > 0) {
    is_valid <- Reduce(`&`, lapply(f_cols, function(col) {
      pred_df[[col]] %in% unique(training_data[[col]])
    }))
    
    pred_df <- pred_df %>%
      dplyr::mutate(dplyr::across(
        .cols = dplyr::all_of(f_cols),
        .fns = ~ factor(.x, levels = unique(training_data[[cur_column()]]))
      ))
    
    pred_df[!is_valid, ] <- NA 
  } else {
    is_valid <- rep(TRUE,nrow(pred_df))
  }
  
  return(list(
    pred_df,
    is_valid
  ))
}
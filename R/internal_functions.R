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
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), if (!is.null(partition) && any(nzchar(partition, keepNA = FALSE))) dplyr::starts_with(partition))
    data <- data.frame(data)
  } else {
    data <- data %>%
      dplyr::select(dplyr::all_of(response), dplyr::all_of(predictors), dplyr::all_of(predictors_f), if (!is.null(partition) && any(nzchar(partition, keepNA = FALSE))) dplyr::starts_with(partition))
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
#wrap_final_list <- function(algo, full_model, variables, response, eval_partial_list, predict_part, part_pred_list, metadata){
wrap_final_list <- function(algo, full_model, variables, response, replica_training_lists, hold_out_evaluation, hold_out_perf, predict_part, metadata){
  # bind predicted evaluation
  eval_partial <- replica_training_lists$eval_partial_list %>%
    dplyr::bind_rows(.id = "replica") %>%
    dplyr::as_tibble()
  
  if (hold_out_evaluation) {
    eval_partial_ho <- replica_training_lists$eval_partial_list_ho %>%
      dplyr::bind_rows(.id = "replica") %>%
      dplyr::as_tibble()
  } else {
    eval_partial_ho <- NULL
  }
  
  # bind predicted partition
  if (predict_part) {
    part_pred <- replica_training_lists$part_pred_list %>%
      dplyr::bind_rows(.id = "replica")
    if (hold_out_evaluation) {
      part_pred_ho <- replica_training_lists$part_pred_list_ho %>%
        dplyr::bind_rows(.id = "replica")
    } else {
      part_pred_ho <- NULL
    }
  } else {
    part_pred <- NULL
    part_pred_ho <- NULL
  }
  
  # Summarize performance
  eval_final <- eval_partial %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(dplyr::across(mae:pdisp, list(
      mean = mean,
      sd = stats::sd
    )), .groups = "drop")
  
  if (hold_out_evaluation) {
    eval_final <- bind_rows(
      eval_final,
      eval_partial_ho %>%
        mutate(model = paste0(model, "_ho")) %>%
        dplyr::group_by(model) %>%
        dplyr::summarise(dplyr::across(mae:pdisp, list(
          mean = mean,
          sd = stats::sd
        )), .groups = "drop")
    )
  }
  
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
    performance_part_hold_out = eval_partial_ho,
    predicted_part = part_pred,
    predicted_part_hold_out = part_pred_ho
  )
  
  # Clean final list
  data_list <- Filter(Negate(is.null), data_list)
  
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
    is_valid <- base::Reduce(`&`, base::lapply(f_cols, function(col) {
      train_levels <- base::levels(training_data[[col]])
      if (base::is.null(train_levels)) train_levels <- base::unique(training_data[[col]])
      
      pred_df[[col]] %in% train_levels
    }))
    
    pred_df[!is_valid, ] <- NA 
    
    pred_df <- pred_df %>%
      dplyr::mutate(dplyr::across(
        .cols = dplyr::all_of(f_cols),
        .fns = function(x) {
          col_name <- dplyr::cur_column()
          
          base::factor(
            x, 
            levels = base::levels(training_data[[col_name]])
          )
        }
      ))
    
  } else {
    is_valid <- rep(TRUE,nrow(pred_df))
  }
  
  return(list(
    pred_df,
    is_valid
  ))
}

check_adapt_holdout_set <- function(
    hold_out_set,
    predictors,
    predictors_f,
    response) {
  if (!is.null(hold_out_set)) {
    hold_out_set$mock_part <- NA
    hold_out_set <- adapt_df(
      data = hold_out_set,
      predictors = predictors,
      predictors_f = predictors_f,
      response = response,
      partition = "mock_part"
    ) %>% select(-mock_part)
  }

  return(hold_out_set)
}

init_training_lists <- function(scopus){
  switch (scopus,
    "replica" = {
      list(
        part_pred_list = list(),
        part_pred_ho_list = list(),
        eval_partial_list = list(),
        eval_partial_ho_list = list()
      )
    },
    "fold" = {
      list(
        eval_partial = list(),
        eval_partial_ho = list(),
        pred_test = list(),
        part_pred = list(),
        part_pred_ho = list()
      )
    }
  )
}

#' fold_perf_register
#'
#' @param model 
#' @param folds 
#' @param j 
#' @param fold_training_lists 
#' @param predict_part 
#' @param hold_out_evaluation 
#' @param pred 
#' @param pred_ho 
#' @param observed 
#' @param observed_ho 
#'
#' @noRd
fold_perf_register <- function(
    model, folds, j,
    fold_training_lists,
    predict_part,
    hold_out_evaluation,
    pred, pred_ho,
    observed, observed_ho) {
  
  fold_training_lists$eval_partial[[j]] <- dplyr::tibble(
    model = model,
    adm_eval(obs = observed, pred = pred)
  )
  
  if (predict_part) {
    fold_training_lists$part_pred[[j]] <- data.frame(partition = folds[j], observed, predicted = pred)
  }
  
  if(hold_out_evaluation){
    fold_training_lists$eval_partial_ho[[j]] <- dplyr::tibble(
      model = model,
      adm_eval(obs = observed_ho, pred = pred_ho)
    )
  }
  
  if(all(predict_part, hold_out_evaluation)){
    fold_training_lists$part_pred_ho[[j]] <- data.frame(partition = folds[j], observed_ho, predicted = pred_ho)
  }
  
  fold_training_lists
}

#' replica_perf_register
#'
#' @param replica_training_lists 
#' @param fold_training_lists 
#' @param folds 
#' @param h 
#' @param predict_part 
#' @param hold_out_evaluation 
#' 
#' @noRd
replica_perf_register <- function(
    replica_training_lists, 
    fold_training_lists,
    folds, 
    h, 
    predict_part, 
    hold_out_evaluation){
  
  names(fold_training_lists$eval_partial) <- 1:length(folds)
  
  fold_training_lists$eval_partial <-
    fold_training_lists$eval_partial[sapply(fold_training_lists$eval_partial, function(x) !is.null(dim(x)))] %>%
    dplyr::bind_rows(., .id = "partition")
  
  replica_training_lists$eval_partial_list[[h]] <- fold_training_lists$eval_partial
  
  if (predict_part) {
    names(fold_training_lists$part_pred) <- 1:length(folds)
    fold_training_lists$part_pred <-
      fold_training_lists$part_pred[sapply(fold_training_lists$part_pred, function(x) !is.null(dim(x)))] %>%
      dplyr::bind_rows(., .id = "partition")
    replica_training_lists$part_pred_list[[h]] <- fold_training_lists$part_pred
  }
  
  if(hold_out_evaluation){
    names(fold_training_lists$eval_partial_ho) <- 1:length(folds)
    
    fold_training_lists$eval_partial_ho <-
      fold_training_lists$eval_partial_ho[sapply(fold_training_lists$eval_partial_ho, function(x) !is.null(dim(x)))] %>%
      dplyr::bind_rows(., .id = "partition")
    
    replica_training_lists$eval_partial_list_ho[[h]] <- fold_training_lists$eval_partial_ho
  }
  
  if(all(hold_out_evaluation, predict_part)){
    names(fold_training_lists$part_pred_ho) <- 1:length(folds)
    fold_training_lists$part_pred_ho <-
      fold_training_lists$part_pred_ho[sapply(fold_training_lists$part_pred_ho, function(x) !is.null(dim(x)))] %>%
      dplyr::bind_rows(., .id = "partition")
    replica_training_lists$part_pred_list_ho[[h]] <- fold_training_lists$part_pred_ho
  }
  
  return(replica_training_lists)
}

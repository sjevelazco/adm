#' Balance dataset at a given presence:absence ratio
#'
#' @param data 
#' @param response 
#' @param absence_ratio 
#'
#' @return
#' @export
#'
#' @examples
balance_dataset <-
  function(data, response, absence_ratio) {
    n_presences <- (data[[response]] > 0) %>% sum
    n_absence = n_presences * absence_ratio
    
    abs_indexes <- which(data[, response] == 0)
    abs_to_remove <- length(abs_indexes) - n_absence
    
    if (abs_to_remove == 0) {
      stop("Dataset already balanced at given ratio.")
    } else if (abs_to_remove < 0) {
      stop("Few absence points to balance at given ratio.")
    }
    
    removed_abs <- sample(abs_indexes, abs_to_remove)
    
    balanced_data <-  data[-removed_abs, ]
    
    return(balanced_data)
  }

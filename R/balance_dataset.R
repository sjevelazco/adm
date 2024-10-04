#' Balance Dataset at a given presence-absence ratio
#'@description The function balances a given dataset based on the specified ratio of presence to absence in the column. 
#' It randomly removes excess of absence in the dataset to achieve the specified ratio.
#' Note that abasence are interpreted all those data with abundance equal to zero.
#' @param data data.frame or tibble. Database that contains a columns with abundance. The dataset to balance.
#' @param response string. The name of the column in data representing the response variable.
#' @param absence_ratio numeric. The desired ratio of presence to absence in the response column.
#' 
#' @return Returns a balanced dataframe with presence:absence ratio in the response column equal to absence_ratio#' 
#' @export
#' 
#' @examples
#' \notrun{
#' require(dplyr)
#' 
#' data("sppabund")
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species three") %>%
#'   dplyr::select(species, ind_ha, x, y)
#' 
#' table(some_sp$ind_ha > 0) 
#' # Note that the dataset is almost balanced
#' # However, as an example, let's assume that we want to reduce 
#' # the number of absences half of the number of presences
#' 
#' some_sp_2 <- balance_dataset(data = some_sp, 
#'                              response = "ind_ha", 
#'                              absence_ratio = 0.5)
#' 
#' table(some_sp$ind_ha > 0)
#' table(some_sp_2$ind_ha > 0)
#' 
#' 
#' 
#' 
#' }
balance_dataset <-
  function(data, response, absence_ratio) {
    n_presences <- (data[[response]] > 0) %>% sum()
    n_absence <- n_presences * absence_ratio

    abs_indexes <- which(data[, response] == 0)
    abs_to_remove <- length(abs_indexes) - n_absence

    if (abs_to_remove == 0) {
      stop("Dataset already balanced at given ratio.")
    } else if (abs_to_remove < 0) {
      stop("Few absence points to balance at given ratio.")
    }

    removed_abs <- sample(abs_indexes, abs_to_remove)

    balanced_data <- data[-removed_abs, ]

    return(balanced_data)
  }

#' Balance database at a given absence-presence ratio
#'
#' @description The function balances a given database based on the specified ratio of absence to presence.
#' It randomly removes excess of absence in the database to achieve the specified ratio.
#' This function interprets as absence all those data with abundance equal to zero.
#'
#' @param data data.frame or tibble. Database that contains a columns with abundance.
#' @param response string. The name of the column in `data` representing the response variable.
#' Note that absence are interpreted all those data with abundance equal to zero.
#' Usage response = "ind_ha"
#' @param absence_ratio numeric. The desired ratio of presence to absence in the
#' response column. E.g., if set to 1 the function will remove absence until have
#' the same number of presence. If set 1.5, the function will remove absence until have
#' 1.5 times the number of presence. Usage absence_ratio = 0.5
#'
#' @return Returns a balanced data.frame or tibble with absence-presence ratio in the response column equal to absence_ratio
#' @export
#'
#' @examples
#' \dontrun{
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
#' some_sp_2 <- balance_dataset(
#'   data = some_sp,
#'   response = "ind_ha",
#'   absence_ratio = 0.5
#' )
#'
#' table(some_sp$ind_ha > 0)
#' table(some_sp_2$ind_ha > 0)
#' }
balance_dataset <-
  function(data, response, absence_ratio) {
    n_presences <- (data[[response]] > 0) %>% sum()
    n_absence <- n_presences * absence_ratio

    abs_indexes <- which(data[, response] == 0)
    abs_to_remove <- length(abs_indexes) - n_absence

    if (abs_to_remove == 0) {
      message("No change in dataset, it is already balanced at given ratio.")
      return(data)
    } else if (abs_to_remove < 0) {
      message("No change in dataset, few absence points to balance at given ratio.")
      return(data)
    }

    removed_abs <- sample(abs_indexes, abs_to_remove)

    balanced_data <- data[-removed_abs, ]

    return(balanced_data)
  }

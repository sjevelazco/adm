#' Performs data transformation on a variable based on the specified method.
#'
#' This function transforms the data in a dataframe or tibble based on the method specified
#' The available methods are "01", "zscore", "log", and "round".
#'
#' @param data A data.frame or tibble containing the data.
#' @param variable A character string specifying the variable (column) to be transformed.
#' @param method A character string specifying the method to be used for transformation. Available methods are "01", "zscore", "log", and "round."
#' \itemize{
#' \item For "01", it scales the variable between 0 and 1 using the formula (x - min(x)) (max(x) - min(x)).
#' \item For "zscore", it standardizes the variable by subtracting the mean and dividing by the standard deviation.
#' \item For "log", it applies natural logarithm transformation to the variable.
#' \item For "round", it rounds the variable's values to the nearest whole numbers.
#' }
#'
#' @return A data.frame or tibble with the transformed variable added as a new column.
#' The new column's name is the original variable name followed by an underscore and method name.
#'
#' @importFrom dplyr pull %>%
#' @export
#'
#' @examples
#' require(dplyr)
#'
#' data("sppabund")
#' # Select data for a single species
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(species, ind_ha, x, y)
#'
#' # Transform abundance data to 0-1
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "01"
#' )
#' some_sp_2
#'
#' # Transform abundance data z-score
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "zscore"
#' )
#' some_sp_2
#'
#' # Transform abundance data log
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "log"
#' )
#' some_sp_2
#'
#' # Transform abundance data log
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "round"
#' )
#' some_sp_2
#'
adm_transform <- function(data, variable, method) {
  x <- data %>% dplyr::pull(variable)

  if (method == "01") {
    # 0 and 1 transformation
    min_ <- min(x, na.rm = TRUE)
    max_ <- max(x, na.rm = TRUE)
    data$.newvar <- ((x - min_) / (max_ - min_))
  } else if (method == "zscore") {
    # zscore transformation
    data$.newvar <- (x - mean(x, na.rm = TRUE)) / (sd(x, na.rm = TRUE))
  } else if (method == "log") {
    # Log transformation
    data$.newvar <- log(x)
  } else if (method == "round") {
    # Round transformation
    data$.newvar <- round(x, digits = 0)
  } else {
    stop("Undefined method")
  }
  names(data)[names(data) == ".newvar"] <-
    paste(variable, method, sep = "_")
  return(data)
}

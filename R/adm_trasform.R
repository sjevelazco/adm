#' Performs data transformation on a variable based on the specified method.
#'
#' This function transforms the data in a tibble or SpatRaster object on the method specified
#' The available methods are "01", "zscore", "log", and "round".
#'
#' @param data A data.frame, tibble, or SpatRaster containing the data.
#' @param variable A character string specifying the variable (column) to be transformed.
#' @param method A character string specifying the method to be used for transformation. Available methods are "01", "zscore", "log", and "round."
#' \itemize{
#' \item For "01", it scales the variable between 0 and 1 using the formula (x - min(x)) (max(x) - min(x)).
#' \item For "zscore", it standardizes the variable by subtracting the mean and dividing by the standard deviation.
#' \item For "log", it applies natural logarithm transformation to the variable.
#' \item For "log1", it sums 1 and then applies natural logarithm transformation to the variable.
#' \item For "round", it rounds the variable's values to the nearest whole numbers.
#' }
#' @param inverse logical. Invert the transformation?
#' @param t_terms vector. c(a,b):
#' \itemize{
#' \item For "01", a = min(x), b = max(x).
#' \item For "zscore", a = mean(x), b = sd(x).
#' \item For "log" and "log1, not needed.
#' \item Can't invert "round" transformations.
#' }
#'
#' @return A data.frame or tibble with the transformed variable added as a new column.
#' The new column's name is the original variable name followed by an underscore and method name.
#'
#' @importFrom dplyr pull %>%
#' @importFrom dplyr pull
#' @importFrom terra global
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' require(terra)
#'
#' # Select data for a single species
#' data("sppabund")
#' some_sp <- sppabund %>%
#'   dplyr::filter(species == "Species one") %>%
#'   dplyr::select(species, ind_ha, x, y)
#'
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)[["bio12"]]
#'
#' # Transform tabular data
#' ## Transform abundance data to 0-1
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "01"
#' )
#' some_sp_2
#'
#' ## Transform abundance data z-score
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "zscore"
#' )
#' some_sp_2
#'
#' ## Transform abundance data log
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "log"
#' )
#' some_sp_2
#'
#' ## Round abundance data
#' some_sp_2 <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "round"
#' )
#' some_sp_2
#'
#' # Tranform raster data
#' ## Transform to 0-1
#' envar_2 <- adm_transform(
#'   data = envar,
#'   variable = "bio12",
#'   method = "01"
#' )
#' envar_2
#'
#' ## Transform z-score
#' envar_2 <- adm_transform(
#'   data = envar,
#'   variable = "bio12",
#'   method = "zscore"
#' )
#' envar_2
#'
#' ## Transform log
#' envar_2 <- adm_transform(
#'   data = envar,
#'   variable = "bio12",
#'   method = "log"
#' )
#' envar_2
#'
#' ## Round
#' envar_2 <- adm_transform(
#'   data = envar,
#'   variable = "bio12",
#'   method = "round"
#' )
#' envar_2
#'
#' # Invert transformation
#' ## Invert 01 tranformation
#' some_sp_transformed <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "01"
#' )
#' some_sp_transformed
#'
#' some_sp_inverted <- adm_transform(
#'   data = some_sp_transformed,
#'   variable = "ind_ha_01",
#'   method = "01",
#'   inverse = TRUE,
#'   t_terms = c(
#'     a = min(some_sp[["ind_ha"]]),
#'     b = max(some_sp[["ind_ha"]])
#'   )
#' )
#' some_sp_inverted
#'
#' ## Invert z-score tranformation
#' some_sp_transformed <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "zscore"
#' )
#' some_sp_transformed
#'
#' some_sp_inverted <- adm_transform(
#'   data = some_sp_transformed,
#'   variable = "ind_ha_zscore",
#'   method = "zscore",
#'   inverse = TRUE,
#'   t_terms = c(
#'     a = mean(some_sp[["ind_ha"]]),
#'     b = sd(some_sp[["ind_ha"]])
#'   )
#' )
#' some_sp_inverted
#'
#' ## Invert log and log1
#' some_sp_transformed <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "log"
#' )
#' some_sp_transformed
#'
#' some_sp_inverted <- adm_transform(
#'   data = some_sp_transformed,
#'   variable = "ind_ha_log",
#'   method = "log",
#'   inverse = TRUE
#' )
#' some_sp_inverted
#'
#' some_sp_transformed <- adm_transform(
#'   data = some_sp,
#'   variable = "ind_ha",
#'   method = "log1"
#' )
#' some_sp_transformed
#'
#' some_sp_inverted <- adm_transform(
#'   data = some_sp_transformed,
#'   variable = "ind_ha_log1",
#'   method = "log1",
#'   inverse = TRUE
#' )
#' some_sp_inverted
#'
#' ## To invert raster, is the same process
#' ## Example:
#' envar_transformed <- adm_transform(
#'   data = envar,
#'   variable = "bio12",
#'   method = "01"
#' )
#' envar_transformed
#'
#' envar_inverted <- adm_transform(
#'   data = envar_transformed,
#'   variable = "bio12_01",
#'   method = "01",
#'   inverse = TRUE,
#'   t_terms = c(
#'     a = terra::global(envar, min, na.rm = T),
#'     b = terra::global(envar, max, na.rm = T)
#'   )
#' )
#' envar_inverted
#' }
adm_transform <- function(data, variable, method, inverse = FALSE, t_terms = NULL) {
  if (inverse) {
    if (class(data)[1] == "SpatRaster") {
      if (method %in% c("zscore", "01") & is.null(t_terms)) {
        stop("Transformation terms are needed to invert zscore and 01 transformations.")
      } else if (method == "round") {
        stop("Round transformation cannot be inverted.")
      } else if (method %in% c("log", "log1")) {
        warning("log and log1 tranformations can't restore original negative values, use carefully.")
      } else if (length(variable) > 1) {
        stop("Can only invert one variable at a time.")
      }

      y <- data[[variable]]

      if (method == "01") {
        # 0 and 1 transformation
        min_ <- t_terms[[1]]
        max_ <- t_terms[[2]]
        x <- (y * (max_ - min_)) + min_
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- t_terms[[1]]
        sd_ <- t_terms[[2]]
        x <- (y * sd_) + mean_
      } else if (method == "log") {
        # Log transformation
        x <- exp(y)
      } else if (method == "log1") {
        # Log + 1 transformation
        x <- exp(y) - 1
      } else {
        stop("Undefined method")
      }

      names(x) <-
        paste(variable, "inverted", sep = "_")

      data <- c(data, x)
    } else if (class(data)[1] %in% c("tbl_df", "df", "data.frame")) {
      if (method %in% c("zscore", "01") & is.null(t_terms)) {
        stop("Transformation terms are needed to invert zscore and 01 transformations.")
      } else if (method == "round") {
        stop("Round transformation cannot be inverted.")
      } else if (method %in% c("log", "log1")) {
        warning("log and log1 tranformations can't restore original negative values, use carefully.")
      } else if (length(variable) > 1) {
        stop("Can only invert one variable at a time.")
      }

      y <- data %>% dplyr::pull(variable)

      if (method == "01") {
        # 0 and 1 transformation
        min_ <- t_terms[[1]]
        max_ <- t_terms[[2]]
        data$.newvar <- (y * (max_ - min_)) + min_
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- t_terms[[1]]
        sd_ <- t_terms[[2]]
        data$.newvar <- (y * sd_) + mean_
      } else if (method == "log") {
        # Log transformation
        data$.newvar <- exp(y)
      } else if (method == "log1") {
        # Log + 1 transformation
        data$.newvar <- exp(y) - 1
      } else {
        stop("Undefined method")
      }
      names(data)[names(data) == ".newvar"] <-
        paste(variable, "inverted", sep = "_")
    }
  } else if (!inverse) {
    if (class(data)[1] == "SpatRaster") {
      x <- data[[variable]]

      if (method == "01") {
        # 0 and 1 transformation
        min_ <- terra::global(x, min, na.rm = TRUE)[[1]]
        max_ <- terra::global(x, max, na.rm = TRUE)[[1]]
        x <- ((x - min_) / (max_ - min_))
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- terra::global(x, mean, na.rm = T)[[1]]
        sd_ <- terra::global(x, sd, na.rm = T)[[1]]
        x <- (x - mean_) / (sd_)
      } else if (method == "log") {
        # Log transformation
        x <- log(x)
      } else if (method == "log1") {
        # Log + 1 transformation
        x <- log(x + 1)
      } else if (method == "round") {
        # Round transformation
        x <- round(x, digits = 0)
      } else {
        stop("Undefined method")
      }

      names(x) <-
        paste(variable, method, sep = "_")

      data <- c(data, x)
    } else if (class(data)[1] %in% c("tbl_df", "df", "data.frame")) {
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
      } else if (method == "log1") {
        # Log + 1 transformation
        data$.newvar <- log(x + 1)
      } else if (method == "round") {
        # Round transformation
        data$.newvar <- round(x, digits = 0)
      } else {
        stop("Undefined method")
      }
      names(data)[names(data) == ".newvar"] <-
        paste(variable, method, sep = "_")
    }
  }

  return(data)
}

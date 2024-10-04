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
adm_transform <- function(data, variable, method, inverse = FALSE, t_terms = NULL) {
  if(inverse){
    if(class(data)[1]=="SpatRaster"){
      if(method %in% c("zscore","01") & is.null(t_terms)){
        stop("Transformation terms are needed to invert zscore and 01 transformations.")
      } else if (method == "round") {
        stop("Round transformation cannot be inverted.")
      } else if (method %in% c("log","log1")) {
        warning("log and log1 tranformations can't restore original negative values, use carefully.")
      } else if (length(variable) > 1){
        stop("Can only invert one variable at a time.")
      }
      
      y <- data[[variable]]
      
      if (method == "01") {
        # 0 and 1 transformation
        min_ <- t_terms[[1]]
        max_ <- t_terms[[2]]
        x <- (y*(max_-min_))+min_
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- t_terms[[1]]
        sd_ <- t_terms[[2]]
        x <- (y*sd_)+mean_
      } else if (method == "log") {
        # Log transformation
        x <- exp(y)
      } else if (method == "log1") {
        # Log + 1 transformation
        x <- exp(y)-1
      } else {
        stop("Undefined method")
      }
      
      names(x) <-
        paste(variable, "inverted", sep = "_")
      
      data <- c(data,x)
    } else if(class(data)[1] %in% c("tbl_df","df","data.frame")) {
      if(method %in% c("zscore","01") & is.null(t_terms)){
        stop("Transformation terms are needed to invert zscore and 01 transformations.")
      } else if (method == "round") {
        stop("Round transformation cannot be inverted.")
      } else if (method %in% c("log","log1")) {
        warning("log and log1 tranformations can't restore original negative values, use carefully.")
      } else if (length(variable) > 1){
        stop("Can only invert one variable at a time.")
      }
      
      y <- data %>% dplyr::pull(variable)
      
      if (method == "01") {
        # 0 and 1 transformation
        min_ <- t_terms[[1]]
        max_ <- t_terms[[2]]
        data$.newvar <- (y*(max_-min_))+min_
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- t_terms[[1]]
        sd_ <- t_terms[[2]]
        data$.newvar <- (y*sd_)+mean_
      } else if (method == "log") {
        # Log transformation
        data$.newvar <- exp(y)
      } else if (method == "log1") {
        # Log + 1 transformation
        data$.newvar <- exp(y)-1
      } else {
        stop("Undefined method")
      }
      names(data)[names(data) == ".newvar"] <-
        paste(variable, "inverted",sep = "_")
    }
    
    
  } else if (!inverse){
    
    if(class(data)[1]=="SpatRaster"){
      x <- data[[variable]]
      
      if (method == "01") {
        # 0 and 1 transformation
        min_ <- terra::global(x,min,na.rm=TRUE)[[1]]
        max_ <- terra::global(x,max,na.rm=TRUE)[[1]]
        x <- ((x - min_) / (max_ - min_))
      } else if (method == "zscore") {
        # zscore transformation
        mean_ <- terra::global(x,mean,na.rm=T)[[1]]
        sd_ <- terra::global(x,sd,na.rm=T)[[1]]
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
      
      data <- c(data,x)
      
    } else if(class(data)[1] %in% c("tbl_df","df","data.frame")) {
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

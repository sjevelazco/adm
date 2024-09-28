#' Extract values from a spatial raster based on x and y coordinates
#'
#' @description The function extracts environmental data at the given x and y coordinates
#' 
#' @param data data.frame or tibble. Database with species abundance with x and y coordinates
#' @param x character. Column name with spatial x coordinates
#' @param y character. Column name with spatial y coordinates
#' @param env_layer SpatRaster. Raster with environmental variables.
#' @param variables character. Vector with the variable names of predictor (environmental) variables
#' Usage variables = c("elevation", "sand", "cfvo"). If no variable is specified, function will
#' return data for all layers. Default NULL
#' @param filter_na logical. If filter_na = TRUE (default), the rows with NA values for any of the
#' environmental variables are removed from the returned tibble.
#'
#' @return
#' A tibble that returns the original data base with additional columns
#' for the extracted environmental variables at each x y location from the SpatRaster object used
#' in 'env_layer'
#'
#' @importFrom dplyr %>% tibble select all_of filter
#' @importFrom stats complete.cases
#' @importFrom terra vect crs extract
#'
#' @export
#'
#' @examples
#' \dontrun{
#' require(terra)
#' 
#' # Datasbase with species abundance and x and y coordinates
#' data("sppabund")
#' 
#' # Raster data with environmental variables
#' envar <- system.file("external/envar.tif", package = "adm")
#' envar <- terra::rast(envar)
#' 
#' # Extract data for a single species
#' some_sp <- sppabund %>%
#'   filter(species == "Species one") %>% dplyr::select(species, ind_ha, x,y)
#' 
#' # Extract environmental data from envar raster for all locations in spp
#' ex_spp <-
#'   adm_extract(
#'     data = some_sp,
#'     x = "x",
#'     y = "y",
#'     env_layer = envar,
#'     variables = NULL,
#'     filter_na = FALSE
#'   )
#' 
#' # Extract environmental for two variables and remove rows with NAs
#' ex_spp2 <-
#'   adm_extract(
#'     data = some_sp,
#'     x = "x",
#'     y = "y",
#'     env_layer = envar,
#'     variables = c("bio1", "elevation"),
#'     filter_na = TRUE
#'   )
#' 
#' ex_spp
#' ex_spp2
#' 
#' 
#' 
#' }
adm_extract <-
  function(data,
           x,
           y,
           env_layer,
           variables = NULL,
           filter_na = TRUE) {

    # Predictor vector when variables.=NULL
    if (is.null(variables)) {
      variables <- names(env_layer)
    }

    # spatial data frame
    sp_data <-
      terra::vect(data,
                  geom = c(x, y),
                  crs = terra::crs(env_layer)
      )

    # extract environmental data at xy locations, if filter_na = FALSE, does not remove rows with NAs
    extract_data <- dplyr::tibble(
      data,
      terra::extract(env_layer[[variables]],
                     sp_data,
                     cells = FALSE
      ) %>%
        dplyr::select({{variables}})
    )


    # removes rows with NAs for any environmental variable
    if (filter_na) {
      complete_vec <- stats::complete.cases(extract_data[, variables])
      if (sum(!complete_vec) > 0) {
        message(sum(!complete_vec), " rows were excluded from database because NAs were found")
        extract_data <- extract_data %>% dplyr::filter(complete_vec)
      }
    }
    return(extract_data)
  }

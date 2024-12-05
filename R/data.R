# data

#' A data set containing species abundance of three species, partition folds, and environmental variables.
#'
#' @format A tibble with 2767 rows and 12 variables:
#' \describe{
#'   \item{species}{species names}
#'   \item{ind_ha}{species abundance expressed as individuals per hectare}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   \item{bio1}{bioclimatic variable related to annual mean temperature}
#'   \item{bio3}{bioclimatic variable related to isothermality}
#'   \item{bio12}{bioclimatic variable related to annual precipitation}
#'   \item{bio15}{bioclimatic variable related to precipitation seasonality}
#'   \item{cfvo}{edaphic variable related to volumetric fraction of coarse fragments}
#'   \item{elevation}{topographic elevation}
#'   \item{sand}{edaphic variable related to soil sand content}
#'   \item{eoc}{ecoregion}
#'   \item{.part1 ... .part3}{repeate k-folds}
#' }
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("sppabund")
#' sppabund
#' }
"sppabund"

#' A data set containing abundance of Cynophalla retusa.
#'
#' @description Cynophalla retusa (Griseb.) Cornejo & Iltis (Capparaceae) is a
#' shrub native to northeastern Argentina, Paraguay, Bolivia,
#' and central Brazil.
#'
#' @format A tibble with 366 rows and 6 variables:
#' \describe{
#'   \item{species}{species names}
#'   \item{ind_ha}{species abundance expressed as individuals per hectare}
#'   \item{x}{longitude of species occurrences}
#'   \item{y}{latitude of species occurrences}
#'   \item{.part}{partition folds}
#' }
#'
#' @examples
#' \dontrun{
#' require(dplyr)
#' data("cretusa_data")
#' cretusa_data
#' }
"cretusa_data"

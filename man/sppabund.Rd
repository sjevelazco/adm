% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sppabund}
\alias{sppabund}
\title{A data set containing species abundance of three species, partition folds, and environmental variables.}
\format{
A tibble with 2767 rows and 12 variables:
\describe{
  \item{species}{species names}
  \item{ind_ha}{species abundance expressed as individuals per hectare}
  \item{x}{longitude of species occurrences}
  \item{y}{latitude of species occurrences}
  \item{bio1}{bioclimatic variable related to annual mean temperature}
  \item{bio3}{bioclimatic variable related to isothermality}
  \item{bio12}{bioclimatic variable related to annual precipitation}
  \item{bio15}{bioclimatic variable related to precipitation seasonality}
  \item{cfvo}{edaphic variable related to volumetric fraction of coarse fragments}
  \item{elevation}{topographic elevation}
  \item{sand}{edaphic variable related to soil sand content}
  \item{eoc}{ecoregion}
  \item{.part1 ... .part3}{repeate k-folds}
}
}
\usage{
sppabund
}
\description{
A data set containing species abundance of three species, partition folds, and environmental variables.
}
\examples{
\dontrun{
require(dplyr)
data("sppabund")
sppabund
}
}
\keyword{datasets}

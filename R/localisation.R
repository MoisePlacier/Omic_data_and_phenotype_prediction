#' Geographic and Population Data for Populus trichocarpa Genotypes
#'
#' This dataset provides geographic coordinates and population identifiers for
#' 200 *Populus trichocarpa* genotypes. Each individual is associated with a
#' specific population label, latitude, longitude, and elevation, along with a
#' numeric class index for use in analyses.
#'
#' @format A data frame with 200 observations of 6 variables:
#' \describe{
#'   \item{Population}{Character — Name of the population or group}
#'   \item{ID}{Character — Unique identifier for each genotype}
#'   \item{Lat}{Numeric — Latitude of the genotype's origin (decimal degrees)}
#'   \item{Lgn}{Numeric — Longitude of the genotype's origin (decimal degrees)}
#'   \item{Elev}{Numeric — Elevation of the origin site (meters above sea level)}
#'   \item{class_index}{Integer — Numeric encoding of the population for analysis}
#' }
#'
#' @source ...
#'
#' @examples
#' # Load the dataset
#' data(localisation)
#'
#' # Summary statistics
#' summary(localisation)
#'
"localisation"

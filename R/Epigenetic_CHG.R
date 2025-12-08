#' CHG DNA Methylation Matrix for Populus trichocarpa Genotypes
#'
#' This dataset contains DNA methylation levels in the CHG context for 199 genotypes
#' of *Populus trichocarpa*, grown in a common garden experiment in Orléans, France.
#' Each value represents the proportion of methylation at a specific CHG site,
#' ranging from 0 (unmethylated) to 1 (fully methylated).
#'
#' @format A numeric matrix with 199 rows and 201,696 columns:
#' \describe{
#'   \item{Rows}{Poplar genotypes (e.g., "1-A02", "1-A04", ..., "1-A25")}
#'   \item{Columns}{Genomic CHG positions (e.g., "chr09_9331153", "chr09_9339147")}
#' }
#'
#' @source ...
#'
#' @examples
#' data(Epigenetic_CHG)
"Epigenetic_CHG"

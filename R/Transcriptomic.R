#' Transcriptomic Expression Matrix for Populus trichocarpa Genotypes
#'
#' This dataset contains normalized gene expression levels for 199 genotypes
#' of *Populus trichocarpa*, grown in a common garden experiment in Orléans, France.
#' Each value represents the expression level of a specific gene in a given genotype,
#' based on RNA-seq data.
#'
#' @format A numeric matrix with 199 rows and 25,819 columns:
#' \describe{
#'   \item{Rows}{Poplar genotypes (e.g., "1-A02", "1-A04", ..., "1-A25")}
#'   \item{Columns}{Gene identifiers from the *Populus trichocarpa* genome v4.1
#'                 (e.g., "Potri_001G000700_2_v4_1", "Potri_001G000800_4_v4_1")}
#' }
#'
#' @source ...
#'
#' @examples
#' data(Transcriptomic)
"Transcriptomic"

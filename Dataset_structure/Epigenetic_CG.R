#' CG DNA Methylation Matrix for Populus trichocarpa Genotypes
#'
#' This dataset contains DNA methylation levels in the CG context for 199 genotypes
#' of *Populus trichocarpa*, grown in a common garden experiment in Orléans, France.
#' Each value represents the proportion of methylation at a specific CG site, ranging
#' from 0 (unmethylated) to 1 (fully methylated).
#'
#' @format A numeric matrix with 199 rows and 313,990 columns:
#' \describe{
#'   \item{Rows}{Poplar genotypes (e.g., "1-A02", "1-A04", ..., "1-A25")}
#'   \item{Columns}{Genomic CG positions (e.g., "chr09_9327121", "chr09_9339362")}
#' }
#'
#' @source ...
#'
#' @examples
#' data(Epigenetic_CG)
"Epigenetic_CG"

#' Kinship Matrix for Populus trichocarpa Genotypes
#'
#' This dataset contains a kinship matrix representing the genetic relatedness
#' among 199 *Populus trichocarpa* genotypes. The matrix is symmetric with
#' dimensions 199 × 199, where each entry quantifies the kinship coefficient
#' between pairs of genotypes.
#'
#' @format A numeric matrix with 199 rows and 199 columns. Row and column names
#' correspond to the unique identifiers of the genotypes.
#'
#' @source ...
#'
#' @examples
#' # Load the kinship matrix
#' data(KinShip_Poplar)
#'
#' # Inspect the matrix dimensions
#' dim(KinShip_Poplar)
#'
#' # View a subset of the kinship matrix
#' KinShip_Poplar[1:5, 1:5]
#'
"KinShip_Poplar"

#' SNP Genotype Matrix for Populus trichocarpa Genotypes
#'
#' This dataset contains SNP genotypes for 199 genotypes of *Populus trichocarpa*,
#' grown in a common garden experiment in Orléans, France.
#' Each value is an integer representing the relationship of the observed nucleotide
#' to the reference genome: 0 = same nucleotide, 1 = complementary base (A<>T, C<>G),
#' 2 = other nucleotide.
#'
#' @format A numeric matrix with 199 rows and 217,000 columns:
#' \describe{
#'   \item{Rows}{Poplar genotypes (e.g., "1-A02", "1-A04", ..., "1-A25")}
#'   \item{Columns}{SNP loci identified by genomic coordinates or marker names
#'   (e.g., "chr01_3051274", "chr08_17359102")}
#' }
#'
#' @source ...
#'
#' @examples
#' # Load the dataset
#' data(Genomic)
#'
#' # Count how many samples are in each group
#' table(Genomic)
#'
"Genomic"

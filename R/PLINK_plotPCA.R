#' PLINK_plotPCA
#' PLINK_plotPCA() generates a PCA plot of the populations from a PLINK bed file
#'
#' @param bed_file  A character string indicating the path of the .bed file
#' @param pop A character string indicating the path of a text file with the
#' population names (identical to the first column in a .fam file)
#'
#' @return pca_plot: ggplot object representing the PCA plot
#' @export
#' @importFrom utils read.table
#' @importFrom BEDMatrix BEDMatrix
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_eig
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' PLINK_pca("data/genotype.bed","data/population.txt")
#' }
#'
PLINK_plotPCA <- function(bed_file,
                      pop){
  # load SNP data using BEDMatrix
  snp_data <- BEDMatrix(bed_file)

  # read family data
  family <- read.table(pop)

  # convert SNP data to a matrix and combine with family data
  snp_matrix <- as.matrix(snp_data)
  combined_data <- cbind(family, snp_matrix)

  # perform PCA
  pca <- PCA(snp_matrix, graph = FALSE)

  # plot individuals with their respective family labels
  pca_plot <- fviz_pca_ind(pca, habillage = combined_data$V1, repel = TRUE)

  return(pca_plot)
}


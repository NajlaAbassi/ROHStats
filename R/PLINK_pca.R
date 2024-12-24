#' PLINK_pca
#' PLINK_pca() generates a PCA plot of the populations from a PLINK bed file
#'
#' @param bed_file  A character string indicating the path of the .bed file
#' @param pop A character string indicating the path of a text file with the
#' population names (identical to the first column in a .fam file)
#'
#' @return pca plot
#' @export
#' @importFrom BEDMatrix BEDMatrix 
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_eig
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' 
PLINK_pca <- function(bed_file, pop){
  # load SNP data using BEDMatrix
  snp_data <- BEDMatrix(bed_file)
  
  # read family data
  family <- read.table(fam_file)
  
  # convert SNP data to a matrix and combine with family data
  snp_matrix <- as.matrix(snp_data)
  combined_data <- cbind(family, snp_matrix)
  
  # perform PCA
  pca <- PCA(snp_matrix, graph = FALSE)
  
  # plot eigenvalues
  fviz_eig(pca)
  
  # plot individuals with their respective family labels
  pca <- fviz_pca_ind(pca, habillage = combined_data$V1, repel = TRUE)
  
  # save plot
  ggsave("figures/pca.png", plot = pca, create.dir = TRUE)
  
  return(pca)
} 


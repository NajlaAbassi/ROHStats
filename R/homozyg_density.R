#' homozyg_density
#'
#' homozyg_density() determines the minimum SNP density (in Mb) for the segment
#' to be considered ROH
#' 
#' @param bim_path a PLINK .bim file
#'
#' @return snp_density
#' @export
#'
#' @examples
#' 
homozyg_density <- function(bim_path) {
  bim_file <- read.table(bim_path, sep = "\t", header = FALSE)
  snp_nb <- nrow(bim_file)
  chromosomes <- unique(bim_file$V1)
  covered_genome_size <- sum(sapply(chromosomes, function(chr) {
    dataset_chr_i <- bim_file[bim_file$V1 == chr, ]
    return(max(dataset_chr_i$V4, na.rm = TRUE, -Inf) - min(dataset_chr_i$V4, na.rm = TRUE, Inf))
  }), na.rm = TRUE)
  density <- snp_nb/covered_genome_size  * 1000
  snp_density <- 1 / density
  return(snp_density)
}



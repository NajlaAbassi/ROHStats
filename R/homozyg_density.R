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
  covered_genome_size = 0
  for (i in 1:22) {
    dataset_chr_i <- bim_file[ bim_file$V1 == i, ]
    size_chr_i <- max(dataset_chr_i$V4, na.rm = FALSE) - min(dataset_chr_i$V4,
                                                             na.rm = FALSE)
    covered_genome_size <- covered_genome_size + size_chr_i
  }
  couverage <- covered_genome_size
  density <- snp_nb/couverage * 1000
  snp_density <- 1 / density
  return(snp_density)
}



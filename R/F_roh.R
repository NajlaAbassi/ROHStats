#' F_roh
#' F_roh() calculates FROH and returns basic statistics on this coefficient
#'
#' @param data_path A character string indicating the path of the .hom.indiv
#' PLINK output
#' @param coverage A character string indicating the path of the .bim file,
#'used to calculate the genome coverage. If the genome size is known
#' it is possible to enter it as a numeric value
#'
#' @return 2 dataframe. The 1st having FROH values and the 2nd having elementary
#' statistics of FROH
#' @export
#'
#' @examples
#' 
F_roh <- function(data_path,coverage) {
  setwd(getwd())
  data <- read.table(data_path, sep = "", header = T)
  id <- as.data.frame(data$IID)
  # If coverage is a file path, read the bim file and compute genome coverage
  if (is.character(coverage) && file.exists(coverage)) {
    bim_file <- read.table(coverage, sep = "\t", header = FALSE)
    snp_nb <- nrow(bim_file)
    ## Compute genome coverage
    covered_genome_size = 0
    for (i in 1:22) {
      dataset_chr_i <- bim_file[bim_file$V1 == i, ]
      size_chr_i <- max(dataset_chr_i$V4, na.rm = FALSE) - min(dataset_chr_i$V4, na.rm = FALSE)
      covered_genome_size <- covered_genome_size + size_chr_i
    }
    coverage <- round(covered_genome_size / 1000000)
  } else {
    # Otherwise, coverage is assumed to be the numeric value `coverage`
    coverage <- as.numeric(coverage)
  }
  # Calculate FROH
  Froh <- data.frame(matrix(NA, nrow = nrow(data), ncol = 1))
  for (i in 1:nrow(data)) {
    Froh[i,1] <- (data$KB[i]/1000) /coverage
  }
  FROH <- cbind(id,Froh)
  colnames(FROH)[1] <- "IID"
  colnames(FROH)[2] <- "Froh"
  FROH$FID <- data$FID
  FROH <- FROH[,c(3,1,2)]
  write.table(FROH,"Froh.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  FROH_values <- read.table("Froh.txt", header = T)
  # Statistics for FROH
  #Mean
  mean_value <- round(mean(FROH$Froh, na.rm = TRUE),2)
  #Standard Deviation
  sd_value <- round(sd(FROH$Froh, na.rm = TRUE),2)
  #Max
  maximum <- max(FROH$Froh)
  #Min
  minimum <- min(FROH$Froh)
  # Summary table
  
  sum_table <- as.data.frame(cbind(maximum,
                                   minimum,
                                   mean_value,
                                   sd_value))
  write.table(sum_table,"Froh_statistics.txt", sep = "\t",col.names = T,
              row.names = F,quote = F)
  summary_table <- read.table("Froh_statistics.txt",header = T)
  return(list(FROH_values = FROH_values, summary_table = summary_table))
}
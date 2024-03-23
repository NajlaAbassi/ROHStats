########## ROH statistical analysis using PLINK outputs ##############
#### Statistics on the number and the length of ROH
####
#' This function performs basic statistics on the number and the length of ROH
#' @param data_path A character string indicating the path of the .hom.indiv PLINK output
#' @return 1 dataframe with elementary statistics
roh_stat <- function(indiv_path) {
  dataset <- read.table(indiv_path, sep = "", header = T)
  ####### total length of ROH per individual##########
  #Maximum
  maximum <- round(max(dataset$KB, na.rm = TRUE)/1000,2)
  #Minimum
  minimum <- round(min(dataset$KB, na.rm = TRUE)/1000,2)
  #Mean
  mean_value <- round(mean(dataset$KB, na.rm = TRUE)/1000,2)
  #Standard Deviation
  sd_value <- round(sd(dataset$KB, na.rm = TRUE)/1000,2)
  # Summary table
  sum_table_1 <- as.data.frame(cbind(maximum,minimum,mean_value,sd_value))
  rownames(sum_table_1)[1] <- "total length"
  write.table(sum_table_1,"total length statistics.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  summary_table_1 <- read.table("total length statistics.txt",header = T)
  ######## number of ROH per individual##########
  #Maximum
  maximum <- round(max(dataset$NSEG, na.rm = TRUE),2)
  #Minimum
  minimum <- round(min(dataset$NSEG, na.rm = TRUE),2)
  #Mean
  mean_value <- round(mean(dataset$NSEG, na.rm = TRUE),2)
  #Standard Deviation
  sd_value <- round(sd(dataset$NSEG, na.rm = TRUE),2)
  # Summary table
  sum_table_2 <- as.data.frame(cbind(maximum,minimum,mean_value,sd_value))
  rownames(sum_table_2)[1] <- "total number"
  write.table(sum_table_2,"mean total length statistics.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  summary_table_2 <- read.table("mean total length statistics.txt",header = T)
  common_col_names <- intersect(names(summary_table_1), names(summary_table_2))
  # Summary table
  summary_table <- merge(summary_table_1, summary_table_2, by = common_col_names, all = T)
  rownames(summary_table)[1] <- "total length"
  rownames(summary_table)[2] <- "total number"
  return(summary_table)
}

###### estimation of genomic inbreeding coefficients FROH and FHOM
######
#' This function calculates FROH and give basic statistics
#' @param data_path A character string indicating the path of the .hom.indiv PLINK output
#' @param bim_path A character string indicating the path of the .bim file, it is used to calculate the genome coverage,
#' if the genome size is known it is possible to enter it as a numeric value
#' @return 2 dataframe; the 1st having FROH values and the 2nd having elementary statistics of FROH
F_roh <- function (data_path,bim_path) {
  setwd(getwd())
  data <- read.table(data_path, sep = "", header = T)
  id <- as.data.frame(data$IID)
  # If bim_path is a file path, read the bim file and compute genome coverage
  if (is.character(bim_path) && file.exists(bim_path)) {
    bim_file <- read.table(bim_path, sep = "\t", header = FALSE)
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
    # Otherwise, bim_path is assumed to be the numeric value `coverage`
    coverage <- as.numeric(bim_path)
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
  write.table(sum_table,"Froh statistics.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  summary_table <- read.table("Froh statistics.txt",header = T)
  return(list(FROH_values = FROH_values, summary_table = summary_table))
}
#####
#' This function calculates FHOM and give basic statistics
#' @param het_path A character string indicating the path of the .het PLINK output
#' @return 2 dataframe; the 1st having FHOM values and the 2nd having elementary statistics of FHOM
F_hom <- function(het_path) {
  dataset<- read.table(het_path, header = T)
  FHOM <- data.frame(dataset$F)
  FHOM$FID <- dataset$FID
  FHOM$IID <- dataset$IID
  FHOM <- FHOM[,c(2,3,1)]
  colnames(FHOM)[3] <- "Fhom"
  write.table(FHOM, 'Fhom.txt', sep = "\t",row.names = FALSE, col.names = T, quote = FALSE)
  FHOM_values <- read.table("Fhom.txt", header = T)
  #Statistics for FHOM
  #Mean
  mean_value <- round(mean(FHOM$Fhom, na.rm = TRUE),2)
  #Standard Deviation
  sd_value <- round(sd(FHOM$Fhom, na.rm = TRUE),2)
  #Max
  maximum <- max(FHOM$Fhom)
  # Min
  minimum <- min(FHOM$Fhom)
  # Summary table
  sum_table <- as.data.frame(cbind(maximum,minimum,mean_value,sd_value))
  write.table(sum_table,"Fhom statistics.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  summary_table <- read.table("Fhom statistics.txt",header = T)
  return(list(FHOM_values = FHOM_values, summary_table = summary_table))
}


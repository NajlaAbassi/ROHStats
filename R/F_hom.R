#' This function calculates FHOM and give basic statistics
#' @param het_path A character string indicating the path of the .het PLINK output
#' @return 2 dataframe; the 1st having FHOM values and the 2nd having elementary statistics of FHOM



#' Title
#'
#' @param het_path A character string indicating the path of the .het PLINK output
#'
#' @return 2 dataframe. The 1st having FHOM values and the 2nd having elementary
#'  statistics of FHOM
#' @export
#'
#' @examples
#' 
F_hom <- function(het_path) {
  dataset <- read.table(het_path, header = T)
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
  write.table(sum_table,"Fhom_statistics.txt", sep = "\t",col.names = T,row.names = F,quote = F)
  summary_table <- read.table("Fhom_statistics.txt",header = T)
  return(list(FHOM_values = FHOM_values, summary_table = summary_table))
}


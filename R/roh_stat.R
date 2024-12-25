#' roh_stat
#' roh_stat() performs basic statistics on the number and the length of ROH
#' @param indiv_path  A character string indicating the path of the .hom.indiv
#' PLINK output
#'
#' @return 1 dataframe with elementary statistics
#' @export
#'
#' @examples
#' 
roh_stat <- function(indiv_path) {
  dataset <- read.table(indiv_path, sep = "", header = T)
  # total length of ROH per individual
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
  write.table(sum_table_1,"total_length_statistics.txt", sep = "\t",
              col.names = T,row.names = F,quote = F)
  summary_table_1 <- read.table("total_length_statistics.txt",header = T)
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
  write.table(sum_table_2,"mean_total_length_statistics.txt", sep = "\t",
              col.names = T,row.names = F,quote = F)
  summary_table_2 <- read.table("mean_total_length_statistics.txt",header = T)
  common_col_names <- intersect(names(summary_table_1), names(summary_table_2))
  # Summary table
  summary_table <- merge(summary_table_1, summary_table_2,
                         by = common_col_names, all = T)
  rownames(summary_table)[1] <- "total length"
  rownames(summary_table)[2] <- "total number"
  return(summary_table)
}
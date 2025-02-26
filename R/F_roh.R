#' F_roh
#' F_roh() calculates FROH and returns basic statistics on this coefficient
#'
#' @param data_path A character string indicating the path of the .hom.indiv
#' PLINK output
#' @param coverage A character string indicating the path of the .bim file, 
#' used to calculate the genome coverage. It is assumed that you have 22 chromosome in this case.
#' If the genome size is known it is possible to enter it as a numeric value
#'
#' @return  @return A list, containing:
#'  \item{FROH_values}{A data frame with FID, IID, and Froh values}
#'  \item{Froh_summary_table}{A data frame with overall Froh summary statistics (Maximum, Minimum, Mean, SD)}
#'  \item{by_population_Froh_stats}{A data frame with per population Froh summary statistics (Maximum, Minimum, Mean, SD)}
#' @export
#' 
#' @import dplyr
#'
#' @examples
#' 
F_roh <- function(data_path,coverage) {
  # read input data
  data <- read.table(data_path, sep = "", header = T)

  # if coverage is a file path, read the bim file and compute genome coverage
  if (is.character(coverage) && file.exists(coverage)) {
    bim_file <- read.table(coverage, sep = "\t", header = FALSE)
    snp_nb <- nrow(bim_file)
    
    ## compute genome coverage
    covered_genome_size <- sum(tapply(bim_file$V4, bim_file$V1, function(x) max(x) - min(x)))
    # convert to megabase
    coverage <- round(covered_genome_size / 1000000)
  } else {
    # Otherwise, coverage is assumed to be the numeric value `coverage`
    coverage <- as.numeric(coverage)
  }
  
  # calculate FROH
  FROH <- data.frame(FID = data$FID, IID = data$IID, Froh = (data$KB / 1000) / coverage)
  
  # Statistics for FROH
  # max
  maximum <- max(FROH$Froh, na.rm = TRUE)
  # min
  minimum <- min(FROH$Froh, na.rm = TRUE)
  # mean
  mean_value <- round(mean(FROH$Froh, na.rm = TRUE),2)
  # standard deviation
  sd_value <- round(sd(FROH$Froh, na.rm = TRUE),2)

  # summary table
  Froh_summary_table <- data.frame(
    Maximum = maximum,
    Minimum = minimum,
    Mean = mean_value,
    SD = sd_value
  )
  
  # compute statistics by population
  by_population_Froh_stats <- FROH %>%
    group_by(FID) %>%
    summarise(
      Maximum = max(Froh, na.rm = TRUE),
      Minimum = min(Froh, na.rm = TRUE),
      Mean = round(mean(Froh, na.rm = TRUE), 2),
      SD = round(sd(Froh, na.rm = TRUE), 2),
      .groups = 'drop'
    ) %>%
    as.data.frame()

  return(list(FROH_values = FROH,
              Froh_summary_table = Froh_summary_table,
              by_population_Froh_stats = by_population_Froh_stats))
}
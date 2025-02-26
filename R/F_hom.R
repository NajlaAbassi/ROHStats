#' F_hom: calculates FHOM and reports basic statistics
#'
#' @param het_path A character string indicating the path of the .het PLINK output
#'
#' @return A list, containing:
#'  \item{FHOM_values}{A data frame with FID, IID, and Fhom values}
#'  \item{Fhom_summary_table}{A data frame with overall Fhom summary statistics (Maximum, Minimum, Mean, SD)}
#'  \item{by_population_Fhom_summary}{A data frame with per population Fhom summary statistics (Maximum, Minimum, Mean, SD)}
#' @export
#' 
#' @import dplyr
#'
#' @examples
#' # simulating a sample PLINK .het dataset
#' sample_data <- data.frame(
#'   FID = c("FAM001", "FAM002", "FAM003"),
#'   IID = c("ID001", "ID002", "ID003"),
#'   O_HOM = c(50, 30, 40),
#'   E_HOM = c(45, 32, 38),
#'   N_NM = c(100, 100, 100),
#'   F = c(0.05, -0.02, 0.08)
#' )
#'
#' # save sample data to a temporary file?
#' temp_file <- tempfile(fileext = ".het")
#' write.table(sample_data, temp_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
#'
#' # run F_hom on the simulated data
#' result <- F_hom(temp_file)
#'
#' # access results
#' result$FHOM_values
#' result$summary_table
#' 
F_hom <- function(het_path) {
  # read input
  dataset <- read.table(het_path, header = T)
  # create Fhom dataframe
  FHOM <- data.frame(FID = dataset$FID, IID = dataset$IID, Fhom = dataset$F)
  # compute stats
  # max
  maximum <- max(FHOM$Fhom, na.rm = TRUE)
  # min
  minimum <- min(FHOM$Fhom, na.rm = TRUE)
  # mean
  mean_value <- round(mean(FHOM$Fhom, na.rm = TRUE),2)
  # standard deviation
  sd_value <- round(sd(FHOM$Fhom, na.rm = TRUE),2)
  
  Fhom_summary_table <- data.frame(
    Maximum = maximum,
    Minimum = minimum,
    Mean = mean_value,
    SD = sd_value
  )
  
  # compute statistics by population
  by_population_Fhom_summary <- FHOM %>%
    group_by(FID) %>%
    summarise(
      Maximum = max(Fhom, na.rm = TRUE),
      Minimum = min(Fhom, na.rm = TRUE),
      Mean = round(mean(Fhom, na.rm = TRUE), 2),
      SD = round(sd(Fhom, na.rm = TRUE), 2),
      .groups = 'drop'
    ) %>%
    as.data.frame()
  
  
  # return output as a list
  return(list(FHOM_values = FHOM,
              Fhom_summary_table = Fhom_summary_table,
              by_population_Fhom_summary = by_population_Fhom_summary))
}


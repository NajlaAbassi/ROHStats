#' relatedness_hm
#'
#' @param data_path A character strings indicating the paths to the .genome file
#' @param ... Additional arguments for pheatmap()
#'
#' @return A heatmap for cryptic relatedness
#' @export
#' @import dplyr
#' @importFrom tidyr pivot_wider
#' @importFrom pheatmap pheatmap
#' 
#'
#' @examples
#' TODO example
relatedness_hm <- function(data_path, ...){
  # read and process the data
  data <- read.table(data_path,header = T)
  
  # convert long-to-wide format
  df_wide <- data %>%
    select(IID1, IID2, PI_HAT) %>%
    tidyr::pivot_wider(names_from = IID2, values_from = PI_HAT)
  
  # convert to matrix with row names
  myM <- as.matrix(df_wide[, -1])
  row.names(myM) <- df_wide$IID1
  
  # convert NAs to 0
  myM[is.na(myM)] <- 0
  
  # plot heatmap
  hm <- pheatmap::pheatmap(myM, ...) 
  
  return(hm)
}

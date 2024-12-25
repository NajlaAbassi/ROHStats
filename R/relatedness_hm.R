#' relatedness_hm
#'
#' @param data_path A character strings indicating the paths to the .genome file
#' @param output_dir A character strings indicating the output directory name.
#' If not indicated it will be "figures/" by default
#'
#' @return A heatmap for cryptic relatedness
#' @export
#' @importFrom reshape2 dcast
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices png
#' @importFrom grid grid.draw
#' 
#'
#' @examples
#' 
relatedness_hm <- function(data_path, output_dir = "figures/"){
  # make sure that the output dir exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # read and process the data
  data <- read.table(data_path,header = T)
  df1 <- data[,c("IID1","IID2","PI_HAT")]
  
  # convert long-to-wide
  x <- reshape2::dcast(df1, IID1 ~ IID2, value.var = "PI_HAT")
  
  # convert to matrix with column and rownames
  myM <- as.matrix(x[ , -1 ])
  row.names(myM) <- x$IID1
  
  # convert all NAs to 0
  myM[ is.na(myM) ] <- 0
  
  # plot heatmap
  hm <- pheatmap::pheatmap(as.data.frame(myM))
  
  # save
  png(filename = file.path(output_dir, "heatmap_pihat.png"),
      width = 10, height = 8, units = "in", res = 300)
  grid::grid.draw(hm$gtable)
  dev.off()
  
  return(hm)
}

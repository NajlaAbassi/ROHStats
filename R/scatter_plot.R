#' scatter_plot
#' scatter_plot() processes multiple input files, associates data with specified
#'  groups, generates individual scatter plots for each group, saves them in
#'  separate files, and creates a combined grid plot for all groups
#'  
#'
#' @param file_paths  A list of character strings indicating the paths to the
#' .hom.indiv files of each group
#' @param groups A character vector of all the groups
#' @param output_dir A character string with the name of the directory where the
#' plot will be saved. It is by default called "figures/"
#'
#' @return individual plots for each group and a combined plot of all the groups
#' @export
#' 
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' 
scatter_plot <- function(file_paths, groups, output_dir = "figures/") {
  # make sure that the output dir exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # read and label
  read_and_label_data <- function(file, group_labels) {
    data <- read.table(file, header = TRUE)
    data <- data[data$FID %in% group_labels, ]
    data$Group <- data$FID
    return(data)
  }
  
  # read and process
  data_list <- lapply(file_paths, read_and_label_data, groups)
  
  # combine all
  total_data <- do.call(rbind, data_list)
  
  # Function to create individual plots
  create_plot <- function(data, title) {
    ggplot(data, aes(x = NSEG, y = KB / 1000)) +
      geom_point() +
      xlim(0, 150) +
      ylim(0, 350) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(y = "Total length of ROH (Mb)", x = "Total number of ROH") +
      ggtitle(title)
  }
  
  # create and save individual plots
  unique_groups <- unique(total_data$Group)
  plots <- list()
  
  for (group in unique_groups) {
    group_data <- total_data[total_data$Group == group, ]
    plot <- create_plot(group_data, group)
    
    # save
    ggsave(
      filename = file.path(output_dir, paste0("plot_length_vs_number_roh_", group, ".png")),
      plot = plot,
      width = 8,
      height = 6
    )
    
    plots[[group]] <- plot
  }
  
  # create a combined plot grid
  num_plots <- length(plots)
  ncol <- ceiling(sqrt(num_plots))
  nrow <- ceiling(num_plots / ncol)
  combined_plot <- grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
  
  # save
  ggsave(
    filename = file.path(output_dir, "combined_plot_length_vs_number_roh.png"),
    plot = combined_plot,
    width = 10,
    height = 8
  )
  
  return(list(individual_plots = plots, combined_plot = combined_plot))
}


#' Generate and Save Grouped Bar Graphs for Each Modified Residue Per Condition
#'
#' @return Grouped bar graphs for each modified residue in the data frame per condition
#' @export
#'
#' @examples generate_grouped_bar_plot_res()
#' @aliases generate_grouped_bar_plot_res
generate_grouped_bar_plot_res <- function() {
  # Prompt the user to select the Excel file
  cat("Select the Excel file containing the quantifiable residue level data: ")
  filepath <- file.choose()

  # Load the data from the selected Excel file
  df_in <- readxl::read_excel(filepath)

  # Auto-detect unique conditions
  unique_conditions <- unique(df_in$Condition)

  # Manually specify blue colors for filling the bars
  condition_colors <- c("#0570b0", "#3690c0", "#74a9cf", "#a6bddb", "#d0d1e6", "#084594", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "#c6dbef", "#08519c", "#3182bd", "#6baed6", "#9ecae1", "#c6dbef")

  # Create a factor variable to represent the sorted order of conditions
  df_in$Condition <- factor(df_in$Condition, levels = unique_conditions)

  # Prompt the user to specify the filename
  cat("Enter the filename for the residue level data (without extension): ")
  filename <- readline()

  # Remove any leading or trailing whitespace
  filename <- trimws(filename)

  # Arrange the dataframe by start
  df_in <- df_in %>%
    arrange(start)

  # Iterate over each protein and make a grouped bar plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # Subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a grouped bar plot of the extent of modification for each residue
    # that maps to this protein, with different conditions represented by color
    fig <- ggplot(temp, aes(x = Res, y = EOM, fill = Condition)) +
      geom_bar(position = position_dodge(width = 0.9), stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), position = position_dodge(width = 0.9), width = 0.25) +
      labs(title = paste(protein, "Residue Level Analysis"),
           x = "Residue",
           y = "Extent of Modification",
           fill = "Condition") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "sans"),
            plot.title = element_text(hjust = 0.5, size = 20, family = "sans", face = "bold"),
            legend.text = element_text(size = 16, family = "sans"),
            legend.title = element_text(size = 18, family = "sans"),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 20))) +
      scale_fill_manual(values = condition_colors)

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(dirname(filepath), paste0(filename, "_ResidueLevelGroupedBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png", width = 10, height = 8, dpi = 300)

    # Print a message to indicate successful saving
    cat("Grouped bar graph for", protein, "saved as", file_out, "\n")
  }
}


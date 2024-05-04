#' Create a Plot for the Extent of Modification at the Peptide Level (Step 19)
#'
#' @param df_in a data frame with quantifiable peptide modifications
#' @param file_output The file path to save the graphs
#' @param filename User defined file name for the saved file
#'
#' @return Bar graphs corresponding to each quantifiable peptide in the input
#' data frame. Must Run save_graphs_pep() to save the graphs in the desired
#' location.
#' @export
#'
#' @examples plot_pep <- generate_eom_plot_pep(df_in = my_data_frame, file_output = "/path/to/output")
#' @aliases generate_eom_plot_pep
#'
generate_eom_plot_pep <- function(df_in, file_output, filename) {
  # Prompt the user to input the Excel file name
  cat("Enter the Excel file name (without extension): ")
  filename <- readline()
  # Remove any leading or trailing whitespace
  filename <- trimws(filename)
  # Get the excel file name from the file path
  filename <- tools::file_path_sans_ext(filename)

  df_in <- df_in %>%
    arrange(start)

  # Create a factor variable to represent the sorted order
  df_in$peptide <- factor(df_in$peptide, levels = df_in$peptide)

  # Iterate over each protein and make an extent of mod plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a bargraph of the extent of modification for each peptide
    # that maps to this protein
    fig <- ggplot(temp, aes(x = peptide, y = EOM)) +
      geom_bar(position = "dodge", stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, linewidth= 1) +
      labs(title = paste(protein," Peptide Level Analysis"),
           x = "Peptide",
           y = "Extent of Modification") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "Arial")) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size = 20,
                                      family = "Arial",
                                      face = "bold")) +
      scale_fill_manual(values = c("grey42"))

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(file_output, paste0(filename, "_PeptideLevelBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png")

    # Print a message to indicate successful saving
    cat("Bar graph for", protein, "saved as", file_out, "\n")
  }
}


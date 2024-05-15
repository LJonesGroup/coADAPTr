#' Create a Plot for the Extent of Modification at the Residue Level (Step 16D)
#'
#' @param df_in a data frame containing the EOM and SD values for each residue
#' @param file_output the directory where the output files will be saved
#' @param filename User defined file name
#'
#' @return a bar graph for each protein with quantifiable residue level
#'  modifications
#' @export
#'
#' @examples  generate_eom_plot_res(df_in = quant_graph_df_res, file_output = file_output,)
#' @aliases generate_eom_plot_res
generate_eom_plot_res <- function(df_in, file_output, filename) {
  # Prompt the user to input the Excel file name
  cat("Enter the Excel file name (without extension): ")
  filename <- readline()
  # Remove any leading or trailing whitespace
  filename <- trimws(filename)
  # Get the excel file name from the file path
  filename <- tools::file_path_sans_ext(filename)

  df_in <- df_in %>%
    arrange(start)

  # Iterate over each protein and make an extent of mod plot for it
  for (protein in unique(df_in$MasterProteinAccessions)) {
    # subset the dataframe for this protein
    temp <- subset(df_in, MasterProteinAccessions == protein)

    # Generate a bargraph of the extent of modification for each residue
    # that maps to this protein
    fig <- ggplot(temp, aes(x = Res, y = EOM)) +
      geom_bar(position = "dodge", stat = "identity") +
      geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, linewidth= 1) +
      labs(title = paste(protein,"Residue Level Analysis"),
           x = "Residue",
           y = "Extent of Modification") +
      theme_classic() +
      theme(text = element_text(size = 18, family = "Arial")) +
      theme(plot.title = element_text(hjust = 0.5,
                                      size = 20,
                                      family = "Arial",
                                      face = "bold")) +
      scale_fill_manual(values = c("grey42"))

    # Create the output directory for bar graphs based on the file output and excel filename
    graph_output_directory <- file.path(file_output, paste0(filename, "_ResidueLevelBarGraphs"))
    dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

    # Generate the full file path for this protein and save the figure
    file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
    ggsave(filename = file_out, plot = fig, device = "png")

    # Print a message to indicate successful saving
    cat("Bar graph for", protein, "saved as", file_out, "\n")
  }
}

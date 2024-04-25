#' generate_eom_plot_res
#'
#' @param df_in a data frame containing the EOM and SD values for each residue
#' @param file_output the directory where the output files will be saved
#' @param excel_filename the name of the excel file that the data was read from
#'
#' @return a bar graph for each protein with quantifiable residue level modification
#' @export
#'
#' @exampleS Bar graphs for each protein with quantifiable residue level modification are generated and saved in a new directory. Must run save_graphs_df_res after to save teh data in a desired location.
generate_eom_plot_res <- function(df_in, file_output, excel_filename) {
  df_in <- df_in %>%
    arrange(start)

  # Grab the protein that is being plotted
  protein <- unique(df_in$MasterProteinAccessions)
  # Generate a bargraph of the extent of modification for each peptide
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = Res, y = EOM)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(aes(ymin = EOM - SD, ymax = EOM + SD), width= .5, linewidth= 1) +
    labs(title = paste(protein," Residue Level Analysis"),
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
  graph_output_directory <- file.path(file_output, paste0(excel_filename, "_ResidueLevelBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

  # Generate the full file path for this protein and save the figure
  file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
  ggsave(filename = file_out, plot = fig, device = "png")

  # Print a message to indicate successful saving
  cat("Bar graph for", protein, "saved as", file_out, "\n")
}

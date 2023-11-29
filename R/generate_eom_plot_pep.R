generate_eom_plot_pep <- function(df_in, file_output, excel_filename) {
  df_in <- df_in %>%
    arrange(start)

  # Create a factor variable to represent the sorted order
  df_in$peptide <- factor(df_in$peptide, levels = df_in$peptide)
  # Grab the protein that is being plotted
  protein <- unique(df_in$MasterProteinAccessions)
  # Generate a bargraph of the extent of modification for each peptide
  # that maps to this protein
  fig <- ggplot(df_in, aes(x = peptide, y = EOM)) +
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
  graph_output_directory <- file.path(file_output, paste0(excel_filename, "_PeptideLevelBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)

  # Generate the full file path for this protein and save the figure
  file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
  ggsave(filename = file_out, plot = fig, device = "png")

  # Print a message to indicate successful saving
  cat("Bar graph for", protein, "saved as", file_out, "\n")
}

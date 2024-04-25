#' generate_eom_plot_pep
#'
#' @param df_in
#' @param file_output
#' @param excel_filename
#'
#' @return Bar graphs corresponding to each quantifiable peptide in the input data frame. Must Run save_graphs_pep() to save the graphs in the desired location.
#' @export
#'
#' @examples plot_pep <- generate_eom_plot_pep(graphing_data = df_in, file_output = file_output, excel_filename = excel_filename)
generate_eom_plot_pep <- function(df_in, file_output, excel_filename) {
  df_in <- df_in %>%
    arrange(start)


  df_in$peptide <- factor(df_in$peptide, levels = df_in$peptide)

  protein <- unique(df_in$MasterProteinAccessions)

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


  graph_output_directory <- file.path(file_output, paste0(excel_filename, "_PeptideLevelBarGraphs"))
  dir.create(graph_output_directory, showWarnings = FALSE, recursive = TRUE)


  file_out <- file.path(graph_output_directory, paste0(protein, ".png"))
  ggsave(filename = file_out, plot = fig, device = "png")


  cat("Bar graph for", protein, "saved as", file_out, "\n")
}

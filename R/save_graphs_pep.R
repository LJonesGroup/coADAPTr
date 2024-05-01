#' save_graphs_pep
#' @param df_in_pep Bar graphs corresponding to peptides that had quantifiable
#' extent of modification
#' @param file_output The file path to save the graphs
#' @param excel_filename The file path to save the excel file
#'
#' @return A series of bar graphs corresponding to peptides that had
#' quantifiable extent of modification
#' @export
#'
#' @examples save_graphs_pep(df_in_pep, file_output, excel_filename)
#' @aliases save_pep_graphs
save_graphs_pep <- function(df, file_output, excel_filename) {
  for (protein in unique(df$MasterProteinAccessions)) {
    temp <- subset(df, MasterProteinAccessions == protein)
    generate_eom_plot_pep(temp, file_output, excel_filename)
  }
}

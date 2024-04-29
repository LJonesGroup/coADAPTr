#' save_graphs_pep
#'
#' @param df_in_pep Bar graphs corresponding to peptides that had quantifiable extent of modification
#'
#' @return A series of bar graphs corresponding to peptides that had quantifiable extent of modification
#' @export
#'
#' @examples
save_graphs_pep <- function(df_in_pep) {
  for (protein in quant_graph_df_pep$MasterProteinAccessions) {
    # subset the dataframe for this protein
    temp <- subset(quant_graph_df_pep, MasterProteinAccessions == protein)
    generate_eom_plot_pep(temp, file_output, excel_filename)
  }
}

#'save_graphs_res
#'
#' @param df
#' @param file_output
#' @param excel_filename
#'
#' @return A series of bar graphs that correspond to the residues with
#' quantifiable modifications
#' @export
#'
#' @examples save_graphs_res(df_in_pep, file_output, excel_filename)
#' @aliases save_res_graphs
save_graphs_res<- function(df, file_output, excel_filename) {
  for (protein in unique(df$MasterProteinAccessions)) {
    temp <- subset(df, MasterProteinAccessions == protein)
    generate_eom_plot_res(temp, file_output, excel_filename)
  }
}

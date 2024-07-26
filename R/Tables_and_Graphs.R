

#' Generate and Save the Tables and Graphs for the LFQ Data
#'
#' @return Tables as excel files and graphs as .PNG files corresponding to the LFQ data.
#' @export
#'
#' @examples Tables_and_Graphs()
#' @aliases Tables_and_Graphs
Tables_and_Graphs<- function() {
  #Create a table of totals
  TotalsTable<<-create_totals_tablelist(graphing_df_pep, graphing_df_res)

  #Save Data Frames as Excel Files
  save_data_frames(file_output, TotalsTable = TotalsTable, Areas_pep = Areas_pep, quant_graph_df_pep = quant_graph_df_pep, Areas_res = Areas_res, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res = graphing_df_res)

  #Generate the grouped bar plots for peptide level data
  generate_grouped_bar_plot_pep()

  #Generate the grouped bar plots for residue level data
  generate_grouped_bar_plot_res()

  #Count the modified peptides and residues per protein
  count_peptides_per_protein()
  count_residue_entries_per_protein()
}


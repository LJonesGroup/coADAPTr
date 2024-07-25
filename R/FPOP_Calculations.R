#' Calculate the Extent of Modification at the Peptide and Residue Level
#'
#' @return A list of dataframes that contain the Extent of Modification (EOM) and SD (Variance) at the peptide and residue level
#' @export
#'
#' @examples FPOP_Calculations()
#' @aliases FPOP_Calculations
FPOP_Calculations<- function() {
  #Calculate the Extent of Modification (EOM) and SD (Variance) at the peptide level
  Areas_pep<<- area_calculations_pep(mod_data_fasta_merged)

  #Merge metadata with numeric graphing data
  graphing_df_pep<<- merge_metadata_pep(Areas_pep, mod_data_fasta_merged)

  #Filter Peptide Level Graphical Data to include data that is acceptable
  quant_graph_df_pep<<- filtered_graphing_df_pep(graphing_df_pep)

  #Calculate the Extent of Modification (EOM) and SD (Variance) at the residue level
  Areas_res<<- area_calculations_res(mod_data_fasta_merged)

  #Subset sequence metadata like residue start/stop for residue level data
  graphing_df_res<<- graphing_data_res(Areas_res, mod_data_fasta_merged)

  #Filter Residue Level Graphical Data to include data that is acceptable
  quant_graph_df_res<<- filtered_graphing_df_res(graphing_df_res)

}


#' Calculate the Extent of Modification at the Peptide and Residue Level
#'
#' @return A list of dataframes that contain the Extent of Modification (EOM) and SD (Variance) at the peptide and residue level
#' @export
#'
#' @examples EOM_Calculations()
#' @aliases EOM_Calculations
EOM_Calculations<- function(EOMbustraction) {

  if(EOMbustraction=="si"){
    #Calculate the Extent of Modification (EOM) and SD (Variance) at the peptide level
    Areas_pep<<- area_calculations_pep(mod_data_fasta_merged, "si")
    write.csv(Areas_pep, file.path(file_output,"Areas_pep.csv"))
  }else{
    #Calculate the Extent of Modification (EOM) and SD (Variance) at the peptide level (no background substraction - not recommended)
    Areas_pep<<- area_calculations_pep(mod_data_fasta_merged, "no")
    write.csv(Areas_pep, file.path(file_output,"Areas_pep_no.csv"))
  }

  #Merge metadata with numeric graphing data
  graphing_df_pep<<- merge_metadata_pep(Areas_pep, mod_data_fasta_merged)

  #Filter Peptide Level Graphical Data to include data that is acceptable
  quant_graph_df_pep<<- filtered_graphing_df_pep(graphing_df_pep)



  if(EOMbustraction=="si"){
    #Calculate the Extent of Modification (EOM) and SD (Variance) at the residue level
    Areas_res<<- area_calculations_res(mod_data_fasta_merged, "si")
  }else{
    #Calculate the Extent of Modification (EOM) and SD (Variance) at the res level (no background substraction - not recommended)
    Areas_res<<- area_calculations_res(mod_data_fasta_merged, "no")
    write.csv(Areas_res, file.path(file_output,"Areas_res_no.csv"))
  }



  #Subset sequence metadata like residue start/stop for residue level data
  graphing_df_res<<- graphing_data_res(Areas_res, mod_data_fasta_merged)

  #Filter Residue Level Graphical Data to include data that is acceptable
  quant_graph_df_res<<- filtered_graphing_df_res(graphing_df_res)

}


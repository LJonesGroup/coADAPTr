for (protein in quant_graph_df_pep$MasterProteinAccessions) {
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_pep, MasterProteinAccessions == protein)
  generate_eom_plot_pep(temp, file_output, excel_filename)
}

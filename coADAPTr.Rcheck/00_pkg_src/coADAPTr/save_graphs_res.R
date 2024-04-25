filtered_graphing_df_pep <- function(df_in_res) {
  for (protein in df_in_res$MasterProteinAccessions) {
  # subset the dataframe for this protein
  temp <- subset(quant_graph_df_res, MasterProteinAccessions == protein)
  generate_eom_plot_res(temp, file_output, excel_filename)
}}

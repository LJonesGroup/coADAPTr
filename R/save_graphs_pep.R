filtered_graphing_df_pep <- function(df_in_pep) {
  for (protein in df_in_pep$MasterProteinAccessions) {
  temp <- subset(quant_graph_df_pep, MasterProteinAccessions == protein)
  generate_eom_plot_pep(temp, file_output, excel_filename)
}}

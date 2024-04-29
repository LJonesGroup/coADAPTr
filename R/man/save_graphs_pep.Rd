\name{save_graphs_pep}
\alias{save_pep_graphs}
\title{Save Extent of Modification (EOM) Plots for modified residues}
\description{
  This function generates extent of modification (EOM) plots for each peptide
  in the input data frame. The plots are saved as PNG files.
}
\usage{
  save_graphs_pep <- function(df_in_pep) {
  for (protein in df_in_pep$MasterProteinAccessions) {
    # subset the dataframe for this protein
    temp <- subset(df_in_pep, MasterProteinAccessions == protein)
    generate_eom_plot_pep(temp, file_output, excel_filename)
  }
}
\arguments{
  \item{df_in_pep}{A data frame containing peptide information for
  plotting.}
  \item{file_output}{The directory path where the graphs will be saved.}
  \item{excel_filename}{The name of the Excel file used to create directories
  for bar graphs.}
}
\details{
  Ensure \code{excel_filename} is defined before using this function.
  It should be generated using \code{tools::file_path_sans_ext(basename(file_path))}.
}
\value{
  Generates extent of modification (EOM) plots for each peptide in the
  input data frame. Saves the graphs as PNG files.
}
\examples{
  excel_filename <- tools::file_path_sans_ext(basename(file_path))
  generate_protein_eom_plots_pep(quant_graph_df_pep, file_output, excel_filename)
}

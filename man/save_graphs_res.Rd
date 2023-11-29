\name{generate_protein_eom_plots_res}
\alias{generate_protein_eom_plots_res}
\title{Generate Extent of Modification (EOM) Plots for Residues}
\description{
  This function generates extent of modification (EOM) plots for each residue in the input data frame.
}
\usage{
  generate_protein_eom_plots_res(quant_graph_df_res, file_output, excel_filename)
}
\arguments{
  \item{quant_graph_df_res}{A data frame containing residue information for plotting.}
  \item{file_output}{The directory path where the graphs will be saved.}
  \item{excel_filename}{The name of the Excel file used to create directories for bar graphs.}
}
\details{
  Ensure \code{excel_filename} is defined before using this function. It should be generated using \code{tools::file_path_sans_ext(basename(file_path))}.
}
\value{
  Generates extent of modification (EOM) plots for each residue in the input data frame. Saves the graphs as PNG files.
}
\examples{
  # Define the excel_filename variable
  excel_filename <- tools::file_path_sans_ext(basename(file_path))

  # Call the function to generate EOM plots for residues
  generate_protein_eom_plots_res(quant_graph_df_res, file_output, excel_filename)
}

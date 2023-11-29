\name{generate_protein_eom_plots}
\alias{generate_protein_eom_plots}
\title{Generate Extent of Modification (EOM) Plots for Proteins}
\description{
  This function generates extent of modification (EOM) plots for each protein in the input data frame.
}
\usage{
  generate_protein_eom_plots(quant_graph_df_pep, file_output, excel_filename)
}
\arguments{
  \item{quant_graph_df_pep}{A data frame containing protein information for plotting.}
  \item{file_output}{The directory path where the graphs will be saved.}
  \item{excel_filename}{The name of the Excel file used to create directories for bar graphs.}
}
\details{
  Ensure \code{excel_filename} is defined before using this function. It is generated using \code{tools::file_path_sans_ext(basename(file_path))}.
}
\value{
  Generates extent of modification (EOM) plots for each protein in the input data frame. Saves the graphs as PNG files.
}
\examples{
  # Define the excel_filename variable
  excel_filename <- tools::file_path_sans_ext(basename(file_path))

  # Call the function to generate EOM plots for proteins
  generate_protein_eom_plots(quant_graph_df_pep, file_output, excel_filename)
}

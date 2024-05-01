\name{save_graphs_pep}
\alias{save_pep_graphs}
\title{Save Extent of Modification (EOM) Plots for modified residues}
\description{
  This function generates extent of modification (EOM) plots for each peptide
  in the input data frame. The plots are saved as PNG files.
}
\usage{
save_graphs_pep <- function(df, file_output, excel_filename) {
  for (protein in unique(df$MasterProteinAccessions)) {
    # subset the dataframe for this protein
    temp <- subset(df, MasterProteinAccessions == protein)
    generate_eom_plot_pep(temp, file_output, excel_filename)
  }
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
 save_graphs_pep(df_in_pep, file_output, excel_filename)
}

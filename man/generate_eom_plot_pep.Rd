\name{generate_peptide_bar_graph}
\alias{generate_peptide_bar_graph}
\title{Generate Peptide Bar Graph}
\description{
  This function generates a bar graph representing the extent of modification for each peptide mapped to a specific protein.
}
\usage{
  generate_peptide_bar_graph(df_in, file_output, excel_filename)
}
\arguments{
  \item{df_in}{A data frame containing information for plotting.}
  \item{file_output}{The directory path where the graphs will be saved.}
  \item{excel_filename}{The name of the Excel file used to create directories for bar graphs.}
}
\value{
  Generates a bar graph representing the extent of modification for each peptide mapped to the protein. Saves the graph as a PNG file.
}
\examples{
  # Call the function to generate peptide bar graphs
  generate_peptide_bar_graph(df_in, file_output, excel_filename)
}

\name{save_data_frames}
\alias{save_df}
\title{Save Data Frames as Separate Excel Files}
\description{
  This function saves multiple data frames as separate Excel files.
}
\usage{
  save_data_frames(file_path, output_directory, ...)
}
\arguments{
  \item{file_path}{The path of the input file.}
  \item{output_directory}{The directory path where the Excel files will be saved.}
  \item{...}{Multiple data frames to be saved.}
}
\details{
  This function extracts the file name from \code{file_path} and removes its extension to generate output file names. It then saves each data frame as a separate Excel file in the \code{output_directory}.
}
\value{
  Saves each provided data frame as a separate Excel file in the specified directory.
}
\examples{
  # Example usage:
  # Define data frames to be saved
  TotalsTable <- data.frame()
  quant_graph_df_pep <- data.frame()
  quant_graph_df_res <- data.frame()
  graphing_df_pep <- data.frame()
  graphing_df_res <- data.frame()

  # Save the data frames as separate Excel files
  save_data_frames(file_path, file_output, TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep, quant_graph_df_res = quant_graph_df_res, graphing_df_pep = graphing_df_pep, graphing_df_res = graphing_df_res)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_data_frames.R
\name{save_data_frames}
\alias{save_data_frames}
\title{Save the Data Frames from your Analysis (Step 16B)}
\usage{
save_data_frames(output_directory, ...)
}
\arguments{
\item{output_directory}{the directory to save the data frames}

\item{...}{data frames to save}
}
\value{
data frames saved in the desired file path
}
\description{
Save the Data Frames from your Analysis (Step 16B)
}
\examples{
save_data_frames(file_path, file_output,
TotalsTable = TotalsTable, quant_graph_df_pep = quant_graph_df_pep,
quant_graph_df_res = quant_graph_df_res)
}

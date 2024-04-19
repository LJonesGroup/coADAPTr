\name{graphing_data_res}
\alias{graphing_data_res}
\title{Prepare Data for Graphing}
\description{
  This function prepares data for graphing by performing various operations, including merging metadata, filtering, sorting, and modifying column values.
}
\usage{
  graphing_data_res(df_in)
}
\arguments{
  \item{df_in}{A data frame containing the input data.}
}
\value{
  Returns a modified data frame suitable for graphing purposes.
}
\examples{
  # Call the function to prepare data for graphing
  graphing_df <- graphing_data_res(Areas_res)
  print(graphing_df)
}

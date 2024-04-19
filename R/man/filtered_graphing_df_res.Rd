\name{filtered_graphing_df_res}
\alias{filtered_graphing_df_res}
\title{Filter Graphing Data}
\description{
  This function filters a data frame for graphing purposes by applying multiple filter conditions.
}
\usage{
  filtered_graphing_df_res(df_in)
}
\arguments{
  \item{df_in}{A data frame containing the input data.}
}
\value{
  Returns a filtered data frame suitable for graphing purposes.
}
\examples{
  # Call the function to filter graphing data
  quant_graph_df_res <- filtered_graphing_df_res(graphing_df_res)
  print(quant_graph_df_res)
}

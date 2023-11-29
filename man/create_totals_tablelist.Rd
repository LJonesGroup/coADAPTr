\name{create_totals_tablelist}
\alias{create_totals_tablelist}
\title{Create Totals Table}
\description{
  This function generates a totals table by counting unique values and quantifiable modifications from input data frames.
}
\usage{
  create_totals_tablelist(df_in, df_res)
}
\arguments{
  \item{df_in}{A data frame containing input data.}
  \item{df_res}{A data frame containing results data.}
}
\value{
  Returns a data frame representing the totals table with counts of unique and quantifiable modifications.
}
\examples{
  # Call the function to create a totals table
  TotalsTable <- create_totals_tablelist(graphing_df_pep, graphing_df_res)
  print(TotalsTable)
}

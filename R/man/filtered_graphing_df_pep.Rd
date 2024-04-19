\name{filtered_graphing_df_pep}
\alias{filter_pep}
\title{
  filter the quantifiable peptide extent of modifications
}
\description{
The function aims to isolate the quantifiable modifications.Specificaly i performs the following:
  It filters rows where the 'EOM' (Extent of Modification) value is greater than zero, indicating a certain level of modification.
  Further, it removes rows where 'EOM' is greater than the standard deviation ('SD'), potentially filtering outliers or specific data points.
  It selectively retains rows where the 'N' column (presumably representing frequency or count) is greater than or equal to 4.
  Finally, it arranges the resulting dataframe based on the 'start' column before returning the refined dataframe.
}
\usage{
  filtered_graphing_df_pep(df)
}
\arguments{
  \item{df}{A data frame consisting of the extent of modification calculations and corresponding colummns. }
}
\examples{
 quant_graph_df_pep <- filtered_graphing_df_pep(graphing_df_pep)
}

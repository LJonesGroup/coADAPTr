\name{grab_seq_metadata_pep}
\alias{grab_pep}
\title{
  grab the quantifiable peptide extent of modifications
}
\description{
  The function aims to isolate the quantifiable modifications.
  It selects specific columns ('MasterProteinAccessions', 'Sequence', 'peptide', 'start') from the input dataframe.
  It removes any duplicated rows based on these selected columns to ensure unique metadata information is retained.
  Finally, it returns a refined dataframe containing unique metadata columns essential for subsequent analysis or operations, discarding any duplicate entries.
}
\usage{
  grab_seq_metadata_pep(df)
}
\arguments{
  \item{df}{A data frame consisting of the extent of modification calculations and corresponding colummns. }
}
\examples{
  graphing_df_pep <- Areas_pep %>%
  left_join(grab_seq_metadata_pep(pd_data_fasta_merged))


graphing_df_pep <- graphing_df_pep[order(graphing_df_pep$start), ]
graphing_df_pep$MasterProteinAccessions <- gsub(".*\\|(.*?)\\|.*", "\\1", graphing_df_pep$MasterProteinAccessions)
}

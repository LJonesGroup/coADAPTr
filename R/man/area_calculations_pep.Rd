\name{area_calculations_pep}
\alias{area_pep}
\title{
  calculate the peptide total areas
}
\description{
  The function aims to prepare and refines the input data for the calculation of the extent of modificaton.
 It standardizes column names by renaming 'Protein Accessions' to 'UniprotID' to match a reference FASTA file.
 It filters out rows containing conserved proteins based on a specific pattern in the 'UniprotID'.
 It cleans 'Modifications' data, removing specific strings like 'Carbamidomethyl' using regular expressions and modifies the data structure to denote the oxidation   status of each entry.
 Lastly, it processes 'Modifications' data to extract positions and types of modifications, creating new columns to store this information.
}
\usage{
  area_calculations_pep(df)
}
\arguments{
  \item{df}{A data frame consisting of the 'raw data' that has been merged with the corresponding annotations from the fasta file }
}
\examples{
  Areas_pep <- area_calculations_pep(pd_data_fasta_merged)
}

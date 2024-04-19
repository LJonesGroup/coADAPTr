\name{parse_fasta}
\alias{parse}
\title{
  parse fasta files
}
\description{
  The function aims to take the input from a FASTA file and modify the column names for storage and usage to map the peptide to the entire protein sequence. Specifically this function performs the following:
  It renames columns within the input dataset, changing 'seq.text' to 'protein_sequence' and 'seq.name' to 'MasterProteinAccessions'.
  Utilizing a regex operation, it extracts and isolates the Uniprot IDs from the 'MasterProteinAccessions' column by pattern matching and substitutions.
  The function primarily standardizes column names and extracts Uniprot IDs from the 'MasterProteinAccessions', preparing the dataset for further analysis or merging with other data sources.
}
\usage{
  parse_fasta(FASTA)
}
\arguments{
  \item{FASTA}{A .FASTA file containing the species data based used to search the MS data. }
}
\examples{
  FASTA <- parse_fasta(FASTA)
}

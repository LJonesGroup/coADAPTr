\name{locate_startend_res}
\alias{locate}
\title{
  locate start and end residues
}
\description{
  The function aims to enrich the 'raw_data' by merging with external data, annotating peptides, and providing details about modifications within the protein sequences.
 Specifically this function performs the following:
 Extract unique entries for 'Master Protein Accessions' and 'Sequence' columns and merge these unique entries with a reference FASTA dataset.
 Use string matching to locate peptide sequences within protein sequences and annotate these peptides with their start and end positions.
 Count modifications, determine modified residues, and combine information to annotate modified amino acids and their positions within the sequences.

}
\usage{
  locate_startend_res(df)
}
\arguments{
  \item{df}{A data frame containing the properly anotated data frame. }
}
\examples{
  pd_data_fasta_merged <- locate_startend_res(pd_data_annotated)
}

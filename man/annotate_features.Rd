\name{annotate_features}
\alias{annotate_features}
\title{
  Annnotate Features
}
\description{
  The function aims to preprocess and modify the input dataset ('raw_data') by standardizing column names, filtering out specific rows, and restructuring content within the 'Modifications' column based on defined patterns and rules. Specifically it does the following:
  It standardizes column names, replacing 'Protein Accessions' with 'UniprotID'.
  Removes rows containing conserved proteins by filtering out Uniprot ID has multple entries.
  Edits the 'Modifications' column, removing instances of 'Carbamidomethyl' tags and standardizing the content format based on specific rules.
  Creates new columns ('MOD' and 'ModPosition') based on the modified 'Modifications' column, classifying and formatting data related to protein modifications.
}
\usage{
  annotate(df)
}
\arguments{
  \item{df}{A dataframe containing sequence search covalent labeling mass spectrometry data that contains the minimum columns Master Protein Accessions, Protein Accession, Sequence, Precursor Abundance, Spectrum File, and Modifications.}
}
\examples{
  pd_data_annotated <- annotate_features(pd_data)
}

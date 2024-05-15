#' Grab Sequence Metadata for Peptide Level Graphing (Step 10)
#' SKIP Using this Function. It is used in merge_metadata_pep function.
#'
#' @param df_in Data frame with original data merged with the FASTA file.
#' locations
#' @return A data frame containing the MasterProteinAccessions, Sequence,
#' peptide, and start columns. ONLY used inside of merge_metadata_pep function.
#' @export
#'
#' @examples new_df <- grab_seq_metadata_pep(raw_data)
#' @aliases grab_seq_metadata_pep
grab_seq_metadata_pep <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence",
                                     "peptide", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}

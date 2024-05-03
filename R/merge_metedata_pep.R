#' Merge Metadata with Peptide Level Graphing Data (Step 11)
#'
#' @param df_in  Data frame with peptide level extent of modification calculations.
#' @param rawdatafastamerged Data frame with raw data and the fasta data merged.
#'
#' @return Data frame with peptide level extent of modification calculations and
#' metadata from the fasta file.
#' @export
#'
#' @examples graphing_df_pep<- merge_metadata_pep(graphing_data_pep, rawdatafastamerged)
#' @aliases merge_metadata_pep
merge_metadata_pep <- function(df_in, rawdatafastamerged) {
  df_out <- left_join(df_in, grab_seq_metadata_pep(rawdatafastamerged))

  df_out <- df_out[order(df_out$start), ]

  df_out$MasterProteinAccessions <- gsub(".*\\|(.*?)\\|.*", "\\1", df_out$MasterProteinAccessions)

  return(df_out)
}

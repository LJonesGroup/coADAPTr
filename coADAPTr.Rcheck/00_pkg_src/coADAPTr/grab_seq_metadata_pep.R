#' grab_seq_metadata_pep
#' @param df_in
#' @return a data frame containing the MasterProteinAccessions, Sequence,
#' peptide, and start columns
#' @export
#'
#' @examples new_df <- grab_seq_metadata_pep(raw_data)
#' @aliases grab_pep
grab_seq_metadata_pep <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence",
                                     "peptide", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}

#' grab_seq_metadata_res
#' @param df_in data frame with original data to match proteins and residue
#' @return A data frame with columns MasterProteinAccessions, Sequence, Res,
#' start
#' @export
#'
#' @examples new_df<-grab_seq_metadata_res(raw_data)
#' @aliases grab_seq_metadata_res
grab_seq_metadata_res <- function(df_in){
  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence", "Res", "start"))
  df_out <- df_out[!duplicated(df_out), ]
  return(df_out)
}

#' graphing_data_res
#' @param df_in data frame with residue level EOM calculations to filter
#' for quantifiable modifications
#' @return a data frame containing quantifiable resiude level modifications
#' @export
#'
#' @examples data_graphing<- graphing_data_res(filtered_graphing_df_res)
#' @aliases graphing_data_res
graphing_data_res <- function(df_in) {
  df_out <- df_in %>%
    left_join(grab_seq_metadata_res(pd_data_fasta_merged)) %>%
    filter(!(is.na(Res) | Res == "")) %>%
    arrange(start) %>%
    mutate(MasterProteinAccessions = gsub(".*\\|(.*?)\\|.*", "\\1", MasterProteinAccessions))

  return(df_out)
}

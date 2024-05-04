#' Filter Data for Graphing the Extent of Modification at the Residue Level (Step 15)
#'
#' @param df_in data frame with residue level EOM calculations to filter
#' for quantifiable modifications
#' @param pd_data_fasta_mergedraw data merged with fasta data
#'
#' @return a data frame containing quantifiable residue level modifications
#' @export
#'
#' @examples data_graphing<- graphing_data_res(filtered_graphing_df_res)
#' @aliases graphing_data_res
graphing_data_res <- function(df_in, pd_data_fasta_merged) {
  df_out <- df_in %>%
    left_join(grab_seq_metadata_res(pd_data_fasta_merged)) %>%
    filter(!(is.na(Res) | Res == "")) %>%
    arrange(start) %>%
    mutate(MasterProteinAccessions = gsub(".*\\|(.*?)\\|.*", "\\1", MasterProteinAccessions))

  return(df_out)
}

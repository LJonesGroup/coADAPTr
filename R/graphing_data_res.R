#' Filter Data for Graphing the Extent of Modification at the Residue Level (Step 14)
#'
#' @param df_in Data frame with residue level EOM calculations and raw data that
#'  was merged with the FASTA file to filter for quantifiable modifications
#' @param pd_data_fasta_mergedraw Raw data merged with FASTA file
#'
#' @return A data frame containing quantifiable residue level modifications
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

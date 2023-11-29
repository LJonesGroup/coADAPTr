graphing_data_res <- function(df_in) {
  df_out <- df_in %>%
    left_join(grab_seq_metadata_res(pd_data_fasta_merged)) %>%
    filter(!(is.na(Res) | Res == "")) %>%
    arrange(start) %>%
    mutate(MasterProteinAccessions = gsub(".*\\|(.*?)\\|.*", "\\1", MasterProteinAccessions))

  return(df_out)
}

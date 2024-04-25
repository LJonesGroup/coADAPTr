grab_seq_metadata_pep <- function(df_in){

  df_out <- subset(df_in, select = c("MasterProteinAccessions", "Sequence",
                                     "peptide", "start"))
  df_out <- df_out[!duplicated(df_out), ]

  return(df_out)
}

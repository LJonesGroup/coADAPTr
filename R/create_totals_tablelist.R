create_totals_tablelist <- function(df_in, df_res) {
  # Initialize an empty data frame
  df_out <- data.frame(
    UniqueProteinMod = NA,
    UniqueSeqMod = NA,
    UniqueResMod = NA,
    QuantifiableModProtein = NA,
    QuantifiableModSeq = NA,
    QuantifiableModRes = NA
  )

  # Count the number of unique MasterProteinAccessions values
  unique_protein_mod <- unique(df_in$MasterProteinAccessions)
  df_out$UniqueProteinMod <- length(unique_protein_mod)

  # Count the number of unique Sequence values
  unique_seq_mod <- unique(df_in$Sequence)
  df_out$UniqueSeqMod <- length(unique_seq_mod)

  # Count the number of unique Residue values
  unique_res_mod <- unique(df_res$Res)
  df_out$UniqueResMod <- length(unique_res_mod)



  # Select rows where EOM is greater than 0 and greater than SD and SD is not NA
  df_filtered <- df_in %>%
    filter(EOM > 0 & EOM > SD & !is.na(SD))
  df_filteredres <- df_res %>%
    filter(EOM > 0 & EOM > SD & SD > 0)

  # Count the number of unique MasterProteinAccessions values in the filtered data frame
  quantifiable_protein_mod <- unique(df_filtered$MasterProteinAccessions)
  df_out$QuantifiableModProtein <- length(quantifiable_protein_mod)


  # Count the number of unique Sequence values in the filtered data frame
  quantifiable_seq_mod <- unique(df_filtered$Sequence)
  df_out$QuantifiableModSeq <- length(quantifiable_seq_mod)


  # Count the number of unique Residue values in the filtered data frame
  quantifiable_res_mod <- unique(df_filteredres$Res)
  df_out$QuantifiableModRes <- length(quantifiable_res_mod)


  # Return the resulting data frame
  return(df_out)
}


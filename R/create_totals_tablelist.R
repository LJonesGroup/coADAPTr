#' Creating a List of Total Peptides and Residues Modified (Step 17)
#'
#' @param df_in Filtered peptide level data frame
#' @param df_res Filtered residue level data frame
#'
#' @return A data frame containing the total number of unique proteins, sequences, and residues modified
#' @export
#'
#' @examples Totals<- create_totals_tablelist(quant_graph_df_pep, quant_graph_df_res)
#' @aliases create_totals_tablelist
create_totals_tablelist <- function(df_in, df_res) {
  df_out <- data.frame(
    UniqueProteinMod = NA,
    UniqueSeqMod = NA,
    UniqueResMod = NA,
    QuantifiableModProtein = NA,
    QuantifiableModSeq = NA,
    QuantifiableModRes = NA
  )


  unique_protein_mod <- unique(df_in$MasterProteinAccessions)
  df_out$UniqueProteinMod <- length(unique_protein_mod)


  unique_seq_mod <- unique(df_in$Sequence)
  df_out$UniqueSeqMod <- length(unique_seq_mod)


  unique_res_mod <- unique(df_res$Res)
  df_out$UniqueResMod <- length(unique_res_mod)


  df_filtered <- df_in %>%
    filter(EOM > 0 & EOM > SD & !is.na(SD))
  df_filteredres <- df_res %>%
    filter(EOM > 0 & EOM > SD & SD > 0)

  quantifiable_protein_mod <- unique(df_filtered$MasterProteinAccessions)
  df_out$QuantifiableModProtein <- length(quantifiable_protein_mod)


  quantifiable_seq_mod <- unique(df_filtered$Sequence)
  df_out$QuantifiableModSeq <- length(quantifiable_seq_mod)


  quantifiable_res_mod <- unique(df_filteredres$Res)
  df_out$QuantifiableModRes <- length(quantifiable_res_mod)


  return(df_out)
}

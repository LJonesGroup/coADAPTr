#' Creating a List of Total Peptides and Residues Modified (Step 16A)
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
    UniqueProteinDet = NA,
    UniqueSeqDet = NA,
    UniqueResDet = NA,
    QuantifiableModProtein = NA,
    QuantifiableModSeq = NA,
    QuantifiableModRes = NA
  )

  # Count unique proteins
  unique_protein_det <- unique(df_in$MasterProteinAccessions)
  df_out$UniqueProteinDet <- length(unique_protein_det)

  # Create a unique identifier for each combination of protein and sequence
  unique_seq_det <- unique(paste(df_in$MasterProteinAccessions, df_in$Sequence, sep = "_"))
  df_out$UniqueSeqDet <- length(unique_seq_det)

  # Create a unique identifier for each combination of protein, sequence, and residue
  unique_res_det <- unique(paste(df_res$MasterProteinAccessions, df_res$Sequence, df_res$Res, sep = "_"))
  df_out$UniqueResDet <- length(unique_res_det)

  # Filtering and quantifying modifications
  df_filtered <- df_in %>%
    dplyr::filter(EOM > 0, EOM > SD, !is.na(SD))
  df_filteredres <- df_res %>%
    dplyr::filter(EOM > 0, EOM > SD, SD > 0)

  # Quantify modifications based on unique proteins
  quantifiable_protein_mod <- unique(df_filtered$MasterProteinAccessions)
  df_out$QuantifiableModProtein <- length(quantifiable_protein_mod)

  # Quantify modifications based on unique sequences (considering the protein)
  quantifiable_seq_mod <- unique(paste(df_filtered$MasterProteinAccessions, df_filtered$Sequence, sep = "_"))
  df_out$QuantifiableModSeq <- length(quantifiable_seq_mod)

  # Quantify modifications based on unique residues (considering both protein and sequence)
  quantifiable_res_mod <- unique(paste(df_filteredres$MasterProteinAccessions, df_filteredres$Sequence, df_filteredres$Res, sep = "_"))
  df_out$QuantifiableModRes <- length(quantifiable_res_mod)

  return(df_out)
}


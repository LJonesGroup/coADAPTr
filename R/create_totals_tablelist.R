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

  # Unique Protein Detected
  df_out$UniqueProteinDet <- length(unique(na.omit(df_in$MasterProteinAccessions)))

  # Unique Sequence Detected (by protein + peptide)
  df_out$UniqueSeqDet <- length(unique(na.omit(paste(df_in$MasterProteinAccessions, df_in$Sequence, sep = "_"))))

  # Unique Residue Detected (by protein + peptide + mod site)
  df_out$UniqueResDet <- length(unique(na.omit(paste(df_res$MasterProteinAccessions, df_res$Sequence, df_res$Res, sep = "_"))))

  # Filter to quantifiable values
  df_filtered <- df_in %>%
    dplyr::filter(!is.na(EOM), !is.na(SD), EOM > 0, SD > 0, EOM > SD)

  df_filteredres <- df_res %>%
    dplyr::filter(!is.na(EOM), !is.na(SD), EOM > 0, SD > 0, EOM > SD)

  # Quantifiable Modifications
  df_out$QuantifiableModProtein <- length(unique(na.omit(df_filtered$MasterProteinAccessions)))
  df_out$QuantifiableModSeq <- length(unique(na.omit(paste(df_filtered$MasterProteinAccessions, df_filtered$Sequence, sep = "_"))))
  df_out$QuantifiableModRes <- length(unique(na.omit(paste(df_filteredres$MasterProteinAccessions, df_filteredres$Sequence, df_filteredres$Res, sep = "_"))))

  return(df_out)
}



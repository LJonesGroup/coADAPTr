#' Filter Data for Creating Extent of Modification Graphs at the Peptide Level (Step 12)
#'
#' @param df_in A data frame with EOM data to filter for quantifiable
#' modifications
#' @return a data frame containing the filtered FPOP data
#' @export
#'
#' @examples graphing_data<- filtered_graphing_df_pep(EOM)
#' @aliases filtered_graphing_df_pep
#'
filtered_graphing_df_pep <- function(df_in) {
  df_out <- df_in[df_in$EOM > 0 & df_in$EOM > df_in$SD & df_in$N > 2, ]
  df_out <- df_out %>%
    arrange(start) %>%
    filter(!is.na(MasterProteinAccessions))

  # Summarize RT to one value per group
  rt_lookup <- mod_data_fasta_merged %>%
    group_by(MasterProteinAccessions, Sequence, Condition) %>%
    summarize(RT = mean(RT, na.rm = TRUE), .groups = "drop")

  # Join summarized RT
  df_out <- left_join(df_out, rt_lookup, by = c("MasterProteinAccessions", "Sequence", "Condition"))

  return(df_out)
}



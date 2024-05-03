#' Filter Data for Creating Extent of Modification Graphs at the Peptide Level (Step 12)
#' @param df_in a data frame with EOM data to filter for quantifiable
#' modifications
#' @return a data frame containing the filtered FPOP data
#' @export
#'
#' @examples graphing_data<- filtered_graphing_df_pep(EOM)
#' @aliases filtered_graphing_df_pep
filtered_graphing_df_pep <- function(df_in) {
  df_out = df_in[df_in$EOM > 0 & df_in$EOM > df_in$SD & df_in$N > 4, ]
  df_out<- df_out %>%
    arrange(start)
  df_out <- df_out[!is.na(df_out$MasterProteinAccessions), ]
  return(df_out)
}


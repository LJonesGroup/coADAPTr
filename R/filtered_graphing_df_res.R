#' Filter Data for Graphing the Extent of Modification at the Residue Level (Step 15)
#'
#' @param df_in A data frame with residue level EOM data to filter for
#' quantifiable residue level modifications
#' @return A data frame containing the filtered data
#' @export
#'
#' @examples graphing_df_res <- filtered_graphing_df_res(graphing_df)
#' @aliases filtered_graphing_df_res
filtered_graphing_df_res <- function(df_in) {
  df_out <- df_in[df_in$EOM > 0 & df_in$EOM > df_in$SD & df_in$N > 3, ]
  df_out <- df_out %>%
    arrange(start) %>%
    filter(!is.na(MasterProteinAccessions))
  return(df_out)
}

#' filtered_graphing_df_pep
#' @param df_in a data frame with EOM data to filter for quantifiable
#' modifications
#' @return a data frame containing the filtered FPOP data
#' @export
#'
#' @examples graphing_data<- filtered_graphing_df_pep(EOM)
#' @aliases filtered_graphing_df_pep
filtered_graphing_df_pep <- function(df_in) {
  df_in <- df_in %>% filter(EOM > 0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>% filter(df_in$N >= 4)  # Filter for the 12th column > 4
  df_in <- df_in %>% arrange(start)
  return(df_in)
}

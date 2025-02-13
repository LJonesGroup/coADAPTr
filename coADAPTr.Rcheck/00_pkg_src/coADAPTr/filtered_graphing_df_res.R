#' filtered_graphing_df_res
#' @param df_in
#' @return A data frame containing the filtered data
#' @export
#'
#' @examples graphing_df_res <- filtered_graphing_df_res(graphing_df)
filtered_graphing_df_res <- function(df_in) {
  df_in <- df_in %>%filter(EOM >0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>% filter(df_in[[12]] > 3)  # Filter for the 12th column > 4
  df_in <- df_in %>% arrange(start)
  return(df_in)
}

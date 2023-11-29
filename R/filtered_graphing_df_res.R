filtered_graphing_df_res <- function(df_in) {
  df_in <- df_in %>%filter(EOM >0)
  df_in <- df_in %>% filter(EOM > SD)
  df_in <- df_in %>% filter(df_in[[12]] > 3)  # Filter for the 12th column > 4
  df_in <- df_in %>% arrange(start)
  return(df_in)
}

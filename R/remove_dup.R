remove_dup <- function(df) {

  identical_cols <- setdiff(names(df), c("Identifying Node", "DeltaScore"))

  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(MultipleID = ifelse(n() > 1, "Yes", "No")) %>%
    ungroup()

  # Group the rows by the columns specified above and find the row with the highest DeltaScore in each group
  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(DMax = ifelse(MultipleID == "Yes" & DeltaScore == max(DeltaScore), "Max", NA_character_),
           Multiple = ifelse(DMax == "Max", "Max", "Min")) %>%
    ungroup()

  # Filter the rows based on the specified criteria (MultipleID == "Yes" and DMax == "Max" or MultipleID == "No")
  df <- df %>%
    filter(MultipleID == "Yes" & DMax == "Max" | MultipleID == "No")

  return(df)
}

#' Remove Duplicate Sequest Node Identifications (Step 5)
#'
#' @param df data frame with original data from PD to remove duplicates
#' created from using the multi-level sequence searching algorithm
#'
#' @return a data frame with no duplicate sequest node identifications
#' @export
#'
#' @examples pd_data <- remove_dup(pd_data)
#' @aliases remove_dup
remove_dup <- function(df) {

  identical_cols <- setdiff(names(df), c("Identifying Node", "DeltaScore"))

  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(MultipleID = ifelse(n() > 1, "Yes", "No")) %>%
    ungroup()

  df <- df %>%
    group_by(across(all_of(identical_cols))) %>%
    mutate(DMax = ifelse(MultipleID == "Yes" & DeltaScore == max(DeltaScore), "Max", NA_character_),
           Multiple = ifelse(DMax == "Max", "Max", "Min")) %>%
    ungroup()

  df <- df %>%
    filter(MultipleID == "Yes" & DMax == "Max" | MultipleID == "No")

  return(df)
}

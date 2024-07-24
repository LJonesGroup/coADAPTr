

#' Interactively Rename Columns Corresponding to TMT Data
#'
#' @param data
#'
#' @return Data frame with the columns renamed by the user
#' @export
#'
#' @examples renamed_data <- rename_columns_interactively(selected_data)
#' @aliases rename_columns_interactively
rename_columns_interactively <- function(data) {
  # Prompt user to select columns for renaming
  selected_columns <- select.list(
    colnames(data),
    multiple = TRUE,
    title = "Select columns to rename:",
    graphics = TRUE
  )

  # Loop through selected columns and rename them
  for (col in selected_columns) {
    new_name <- readline(prompt = paste("Enter new name for column", col, "in the format Condition:SampleType(L or NL): "))
    # Rename the column
    colnames(data)[colnames(data) == col] <- new_name
  }

  return(data)
}



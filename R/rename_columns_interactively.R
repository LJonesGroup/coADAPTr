

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
    title = "Select columns to rename (Hold Ctrl or Cmd to select multiple):",
    graphics = TRUE
  )

  # Exit if no columns are selected
  if (length(selected_columns) == 0) {
    cat("No columns selected. Exiting without renaming.\n")
    return(data)
  }

  # Loop through selected columns and rename them
  for (col in selected_columns) {
    prompt_msg <- paste0(
      "Enter new name for column '", col, "' in the format ",
      "Condition:SampleType (e.g., Drug:NL1 or Drug Free:L3).\n",
      "For multiplexed experiments, 'NL' (Non Irradiated) samples include both baseline controls (no reagents + no irradiation) ",
      "and no-laser samples with reagents:\n"
    )
    new_name <- readline(prompt = prompt_msg)

    # Trim whitespace and apply the new name
    new_name <- trimws(new_name)

    if (nzchar(new_name)) {
      colnames(data)[colnames(data) == col] <- new_name
    } else {
      cat("No new name entered for column '", col, "'. Skipping rename.\n")
    }
  }

  return(data)
}




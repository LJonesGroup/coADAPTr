
#' Calculate the Sum of Abundances and True Abundances for TMT Data
#'
#' @param data
#'
#' @return A data frame with the sum of abundances and true abundances for TMT data
#' @export
#'
#' @examples TMT_data <- sum_and_calculate_abundancesTMT(raw_data)
#' @aliases sum_and_calculate_abundancesTMT
sum_and_calculate_abundancesTMT <- function(data) {
  # Function to prompt for column selection within the main function
  select_tmt_abundance_columns <- function() {
    col_names <- colnames(data)
    selected_cols <- select.list(
      col_names,
      multiple = TRUE,
      title = "Select columns containing TMT abundances (Hold Ctrl to select multiple):",
      graphics = TRUE
    )
    if (length(selected_cols) == 0) {
      stop("No columns selected, exiting the function.")
    }
    return(selected_cols)
  }

  # Prompt user to select TMT abundance columns
  cat("Please select the columns containing TMT abundances. You can select multiple columns.\n")
  abundance_cols <- select_tmt_abundance_columns()

  # Sum numerical content of selected abundance columns for each row
  data$AbundanceTotals <- rowSums(data[, abundance_cols, drop = FALSE], na.rm = TRUE)

  # Loop through each selected abundance column to calculate true abundances
  for (col in abundance_cols) {
    # Create a new column name based on the original one
    new_col_name <- paste("TRUE Abundance", col, sep = ": ")
    data[[new_col_name]] <- data[[col]] * data$Intensity / data$AbundanceTotals
  }

  return(data)
}

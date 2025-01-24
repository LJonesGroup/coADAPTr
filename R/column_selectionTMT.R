

#' Select the Columns Necessary for Extent of Modification Calculation from TMT Data
#'
#' @param df
#'
#' @return A data frame with selected columns from the input TMT data frame
#' @export
#'
#' @examples selected_data <- column_selectionTMT(TMT_data)
#' @aliases column_selectionTMT
column_selectionTMT <- function(df) {
  refined_data <- data.frame(matrix(ncol = 0, nrow = nrow(df)))  # Initialize an empty data frame with the same number of rows as df
  col_names <- colnames(df)

  # Function to prompt for column selection, allowing multiple selections
  select_column <- function(prompt, allow_multiple = FALSE) {
    cat(prompt, "\n")
    selected_cols <- select.list(
      col_names,
      multiple = allow_multiple,
      title = prompt,
      graphics = TRUE
    )
    if (length(selected_cols) == 0) {
      stop("No column selected, exiting the function.")
    }
    return(selected_cols)
  }

  # Prompt and select columns
  seq_col <- select_column("Please select the column containing the peptide sequences (unannotated):")
  acc_col <- select_column("Please select the column containing the Uniprot ID or master protein accessions:")
  mod_col <- select_column("Please select the column containing the protein modifications:")
  rt_col <- select_column("Please select the column containing the sequence Retention Time:")
  prec_cols <- select_column("Please select the column(s) containing the adjusted precursor abundances/intensities (Hold Ctrl to select multiple):", allow_multiple = TRUE)

  # Add selected columns to refined_data
  refined_data <- cbind(refined_data, df[[seq_col]], df[[acc_col]], df[[mod_col]], df[[rt_col]], df[, prec_cols])

  # Rename columns for clarity if needed
  # Adjust this part to dynamically handle multiple precursor columns
  new_col_names <- c("Sequence", "Master Protein Accessions", "Modifications", "RT")
  new_col_names <- c(new_col_names, prec_cols)  # Append precursor column names directly
  colnames(refined_data) <- new_col_names

  return(refined_data)
}


#' Column Selection for LFQ Data (Non TMT Labeled)
#'
#' @param df
#'
#' @return A data frame with the selected columns.
#' @export
#'
#' @examples Selected_data <- column_selectionLFQ(raw_data)
#' @aliases column_selectionLFQ
column_selectionLFQ <- function(df) {
  refined_data <- data.frame(matrix(ncol = 0, nrow = nrow(df)))  # Initialize an empty data frame with the same number of rows as df
  col_names <- colnames(df)

  # Function to prompt for column selection
  select_column <- function(prompt) {
    cat(prompt, "\n")
    selected_col <- select.list(
      col_names,
      multiple = FALSE,
      title = prompt,
      graphics = TRUE
    )
    if (!is.null(selected_col) && selected_col != "") {
      return(selected_col)
    } else {
      stop("No column selected, exiting the function.")
    }
  }

  # Prompt and select columns
  seq_col <- select_column("Please select the column containing the peptide sequences (unannotated):")
  acc_col <- select_column("Please select the column containing the Uniprot ID or master protein accessions:")
  mod_col <- select_column("Please select the column containing the protein modifications:")
  pre_col <- select_column("Please select the column containing precursor abundance/intensities:")
  spe_col <- select_column("Please select the column containing the spectrum file IDs:")
  RT <- select_column("Please select the column containing the retention time:")

  # Add selected columns to refined_data
  refined_data <- cbind(refined_data, df[[seq_col]], df[[acc_col]], df[[mod_col]], df[[pre_col]], df[[spe_col]], df[[RT]])

  # Rename columns
  colnames(refined_data) <- c("Sequence", "Master Protein Accessions", "Modifications", "Precursor Abundance", "Spectrum File", "RT")

  return(refined_data)
}

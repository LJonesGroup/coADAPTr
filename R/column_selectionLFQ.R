
#' Column Selection for LFQ Data (Non TMT Labeled)
#'
#' @param df
#'
#' @return A data frame with the selected columns.
#' @export
#'
#' @examples Selected_data <- column_selectionLFQ(raw_data)
#' @aliases column_selectionLFQ
column_selectionLFQ <- function(df, filetype) {


  if(filetype == "psm"){

    # the peptide sequences (unannotated)
    # So why not df...it does not work, plus in LFQ_prerp raw data is used as a global environment variable already
    seq_col <- raw_data$'Peptide'
    # the Uniprot ID or master protein accessions
    acc_col <- raw_data$'Protein ID'
    # Modifications
    mod_col <- raw_data$'FPOP Modifications'
    # Precursor abundance/intensities
    pre_col <- raw_data$'Intensity'
    # The spectrum file IDs
    spe_col <- raw_data$'Spectrum'

    refined_data <- data.frame(matrix(ncol = 0, nrow = nrow(raw_data)))  # Initialize an empty data frame with the same number of rows as df
    col_names <- colnames(raw_data)

    refined_data$'Peptide' <- seq_col
    refined_data$'Protein ID' <- acc_col
    refined_data$'FPOP Modifications' <- mod_col
    refined_data$'Intensity' <- pre_col
    refined_data$'Spectrum' <- spe_col



  }else if(filetype=="report"){

    # the peptide sequences (unannotated)
    # So why not df...it does not work, plus in LFQ_prerp raw data is used as a global environment variable already
    seq_col <- raw_data$`Stripped Sequence`
    # the Uniprot ID or master protein accessions
    acc_col <- raw_data$'Protein Ids'
    # Modifications
    mod_col <- raw_data$'FPOP Modifications'
    # Precursor abundance/intensities
    pre_col <- raw_data$`Ms1 Area`
    # The spectrum file IDs
    spe_col <- raw_data$'Run'

    refined_data <- data.frame(matrix(ncol = 0, nrow = nrow(raw_data)))  # Initialize an empty data frame with the same number of rows as df
    col_names <- colnames(raw_data)

    refined_data$'Peptide' <- seq_col
    refined_data$'Protein ID' <- acc_col
    refined_data$'FPOP Modifications' <- mod_col
    refined_data$'Intensity' <- pre_col
    refined_data$'Spectrum' <- spe_col
  }else{
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

    # Add selected columns to refined_data
    refined_data <- cbind(refined_data, df[[seq_col]], df[[acc_col]], df[[mod_col]], df[[pre_col]], df[[spe_col]])



  }

  # Rename columns
  colnames(refined_data) <- c("Sequence", "Master Protein Accessions", "Modifications", "Precursor Abundance", "Spectrum File")

  #print(refined_data)

  return(refined_data)
}

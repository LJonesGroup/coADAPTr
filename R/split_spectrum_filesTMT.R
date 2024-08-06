
#' Extract Spectrum File Identifiers for TMT Data
#'
#' @param df_in
#'
#' @return A data frame containing the split spectrum files.
#' @export
#'
#' @examples mod_data_fasta_merged <- split_spectrum_filesTMT(mod_data_fasta_merged)
#' @aliases split_spectrum_filesTMT
split_spectrum_filesTMT <- function(df_in) {
  # Check if the "Spectrum File" column exists
  if ("Spectrum File" %in% colnames(df_in)) {
    # Initialize the "Condition" column
    df_in$Condition <- NA

    # Loop through each row and split the "Spectrum File" contents
    df_in <- df_in %>%
      mutate(
        Condition = sapply(strsplit(as.character(`Spectrum File`), ":"), `[`, 1),
        `Spectrum File` = sapply(strsplit(as.character(`Spectrum File`), ":"), `[`, 2)
      )
  } else {
    cat("'Spectrum File' column not found in the data frame.\n")
  }

  return(df_in)
}

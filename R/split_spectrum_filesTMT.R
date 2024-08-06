
#' Extract Spectrum File Identifiers for TMT Data
#'
#' @param df_in
#'
#' @return A data frame with the split spectrum files.
#' @export
#'
#' @examples mod_data_fasta_merged <- split_spectrum_filesTMT(mod_data_fasta_merged)
#' @aliases split_spectrum_filesTMT
split_spectrum_filesTMT <- function(df_in) {
  # Check if the "Spectrum File" column exists
  if ("Spectrum File" %in% colnames(df_in)) {
    # Extract the part before the first '.' to find unique files
    df_in$FileIdentifier <- sapply(strsplit(as.character(df_in$`Spectrum File`), "\\."), `[`, 1)
    unique_files <- unique(df_in$FileIdentifier)
    cat("Unique File Identifiers detected:\n")
    print(unique_files)

    # Initialize a data frame to hold the mappings
    mappings <- data.frame(Original = character(), NewName = character(), Condition = character(), stringsAsFactors = FALSE)

    # Loop through each unique file identifier
    for (file in unique_files) {
      response <- readline(prompt = paste("Enter new name for '", file, "' in the format 'Condition:SampleType(L or NL)': ", sep=""))
      parts <- strsplit(response, ":")[[1]]
      if (length(parts) == 2) {
        # Append the mappings
        mappings <- rbind(mappings, data.frame(Original = file, NewName = parts[2], Condition = parts[1]))
      } else {
        cat("Invalid input. Skipping '", file, "'\n")
      }
    }

    # Rename and assign based on mappings
    if (nrow(mappings) > 0) {
      # Map FileIdentifier to NewName and Condition
      for (i in 1:nrow(mappings)) {
        df_in$`Spectrum File`[df_in$FileIdentifier == mappings$Original[i]] <- mappings$NewName[i]
        df_in$Condition[df_in$FileIdentifier == mappings$Original[i]] <- mappings$Condition[i]
      }
      # Optionally remove the temporary identifier column
      df_in$FileIdentifier <- NULL
    }
  } else {
    cat("'Spectrum File' column not found in the data frame.\n")
  }

  return(df_in)
}

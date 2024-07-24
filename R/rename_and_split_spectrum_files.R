
#' Rename and Extract Spectrum File Identifiers
#'
#' @param df_in
#'
#' @return A data frame with the renamed and split spectrum files.
#' @export
#'
#' @examples modified_data <- rename_and_split_spectrum_files(Selected_data)
#' @aliases rename_and_split_spectrum_files
rename_and_split_spectrum_files <- function(df_in) {
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

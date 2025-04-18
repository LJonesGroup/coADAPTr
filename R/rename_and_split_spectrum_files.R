
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
    mappings <- data.frame(
      Original = character(),
      NewName = character(),
      Condition = character(),
      stringsAsFactors = FALSE
    )

    # Loop through each unique file identifier
    for (file in unique_files) {
      response <- readline(
        prompt = paste0(
          "Enter new name for '", file, "'\n",
          "Format: Condition:SampleType (e.g., Drug:NL1 or Drug Free:L3)\n",
          "For multiplexed experiments:\n",
          "- 'NL' (Non Irradiated) includes:\n",
          "  • Baseline control (no reagents AND no irradiation)\n",
          "  • No-laser samples with reagents\n",
          "Your input: "
        )
      )

      parts <- strsplit(trimws(response), ":")[[1]]
      if (length(parts) == 2) {
        # Append the mappings
        mappings <- rbind(mappings, data.frame(
          Original = file,
          NewName = trimws(parts[2]),
          Condition = trimws(parts[1]),
          stringsAsFactors = FALSE
        ))
      } else {
        cat("Invalid input. Skipping '", file, "'\n")
      }
    }

    # Make sure 'Condition' column exists
    if (!"Condition" %in% colnames(df_in)) {
      df_in$Condition <- NA
    }

    # Rename and assign based on mappings
    for (i in seq_len(nrow(mappings))) {
      df_in$`Spectrum File`[df_in$FileIdentifier == mappings$Original[i]] <- mappings$NewName[i]
      df_in$Condition[df_in$FileIdentifier == mappings$Original[i]] <- mappings$Condition[i]
    }

    # Remove the temporary identifier column
    df_in$FileIdentifier <- NULL
  } else {
    cat("'Spectrum File' column not found in the data frame.\n")
  }

  return(df_in)
}


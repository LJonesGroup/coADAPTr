#'Perform Annotation of Raw Data (Step 6)
#'
#' @param raw_data Raw data from proteome discoverer that has the Sample and Control
#' files indicated
#' @return A modified data frame with appropriate modification to the raw data
#' @export
#'
#' @examples newdf<- annotate_features(raw_data)
#' @aliases annotate_features
annotate_features <- function(raw_data) {
  raw_data <- raw_data %>%
    mutate(`Master Protein Accessions` = sapply(strsplit(`Master Protein Accessions`, ";"), `[`, 1),
           UniprotID = `Master Protein Accessions`)  # Adding this line to create the UniprotID column

  # Rename Modifications column to Mods
  raw_data <- raw_data %>%
    rename(Mods = Modifications)

  # Remove specific modifications and clean up semicolons
  raw_data <- raw_data %>%
    mutate(Modifications = Mods) %>%
    mutate(Modifications = gsub(";[A-Z]\\d+\\(Carbamidomethyl\\)", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\);", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\)", "", Modifications)) %>%
    mutate(Modifications = gsub(";[A-Z]\\d+\\(TMT[^)]*\\)", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(TMT[^)]*\\);", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(TMT[^)]*\\)", "", Modifications))

  # Cleanup: Remove leading, trailing, and multiple consecutive semicolons
  raw_data <- raw_data %>%
    mutate(Modifications = gsub("\\s*;\\s*", ";", Modifications)) %>%  # Remove spaces around semicolons
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%  # Replace multiple semicolons with a single one
    mutate(Modifications = gsub("^;|;$", "", Modifications))  # Remove leading and trailing semicolons

  # Determine if the sequence is oxidized or unoxidized
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(is.na(Modifications) | Modifications == "", "Unoxidized",
                        ifelse(grepl("Oxidation", Modifications), "Oxidized", "Unoxidized")))

  # Create ModPositionL and ModPositionN columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = sub("^([A-Z]*).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = gsub(".*?([0-9]+).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = ifelse(ModPositionN == Modifications, NA, ModPositionN)) %>%
    mutate(ModPosition = ifelse(is.na(ModPositionL) | ModPositionL == "", NA,
                                paste(ModPositionL, ModPositionN, sep = "")))

  # Extract letters before ":" in the "Spectrum File" column and add to "Condition" column
  raw_data <- raw_data %>%
    mutate(Condition = sub("^(.*):.*", "\\1", `Spectrum File`))

  return(raw_data)
}


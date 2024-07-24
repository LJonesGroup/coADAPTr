#'Perform Annotation of Raw Data (Step 6)
#'
#' @param raw_data Raw data from proteome discoverer that has the Sample and Control
#' files indicated
#' @return A modified data frame with appropriate modification to the raw data
#' @export
#'
#' @examples modified_data_annotated <- annotate_features(modified_data)
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
    mutate(Modifications = gsub("[A-Z]\\d+\\(TMT[^)]*\\)", "", Modifications)) %>%
    mutate(Modifications = gsub("N-Term\\(Prot\\)\\(TMTpro\\);", "", Modifications)) %>%
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%
    mutate(Modifications = gsub("^;|;$", "", Modifications))

  # Cleanup: Remove leading, trailing, and multiple consecutive semicolons
  raw_data <- raw_data %>%
    mutate(Modifications = gsub("\\s*;\\s*", ";", Modifications)) %>%
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%
    mutate(Modifications = gsub("^;|;$", "", Modifications))

  # Determine if the sequence is oxidized or unoxidized
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(is.na(Modifications) | Modifications == "" | Modifications == "NA" | Modifications == " ", "Unoxidized", "Oxidized"))

  # Create ModPositionL and ModPositionN columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = ifelse(nchar(Modifications) > 0, gsub(".*?([A-Za-z]).*", "\\1", Modifications), NA)) %>%
    mutate(ModPositionN = ifelse(nchar(Modifications) > 0, gsub(".*?(\\d+).*", "\\1", Modifications), NA))

  return(raw_data)
}



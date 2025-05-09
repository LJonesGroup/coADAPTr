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
  library(dplyr)
  library(stringr)

  # Clean up Uniprot ID
  raw_data <- raw_data %>%
    mutate(`Master Protein Accessions` = sapply(strsplit(`Master Protein Accessions`, ";"), `[`, 1),
           UniprotID = `Master Protein Accessions`)

  # Rename and normalize the modifications column
  raw_data <- raw_data %>%
    rename(Mods = Modifications) %>%
    mutate(Modifications = Mods) %>%
    # Replace commas with semicolons for uniform parsing
    mutate(Modifications = str_replace_all(Modifications, ",", ";")) %>%
    # Remove Carbamidomethyl on C based on mass (57.*)
    mutate(Modifications = str_remove_all(Modifications, "\\d*C\\(57(?:\\.\\d+)?\\)")) %>%
    # Remove TMT-related modifications
    mutate(Modifications = str_remove_all(Modifications, "[A-Z]\\d+\\(TMT[^)]*\\)")) %>%
    # Remove N-Term(TMT) or N-term modifications
    mutate(Modifications = str_remove_all(Modifications, "N[-|_]term\\([^)]*\\)")) %>%
    # Clean extra separators
    mutate(Modifications = str_replace_all(Modifications, "\\s*;\\s*", ";")) %>%
    mutate(Modifications = str_replace_all(Modifications, ";{2,}", ";")) %>%
    mutate(Modifications = str_remove_all(Modifications, "^;|;$"))

  # Flag as oxidized/unoxidized
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(Modifications == "" | is.na(Modifications) | Modifications == "NA",
                        "Unoxidized", "Oxidized"))

  # Extract first biologically meaningful mod (e.g., 9H, 8E)
  raw_data <- raw_data %>%
    mutate(first_mod_entry = str_extract(Modifications, "\\d+[A-Z]\\([^)]*\\)"),
           first_mod = str_extract(first_mod_entry, "[A-Z]"),
           first_pos = str_extract(first_mod_entry, "\\d+"))

  # Final columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = first_mod,
           ModPositionN = as.numeric(first_pos))

  return(raw_data)
}

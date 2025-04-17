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
    # Remove common unwanted tags (Carbamidomethyl, TMT tags, etc.)
    mutate(Modifications = str_remove_all(Modifications, "[A-Z]\\d+\\((Carbamidomethyl|TMT[^)]*)\\)")) %>%
    mutate(Modifications = str_remove_all(Modifications, "N-Term\\(Prot\\)\\(TMTpro\\)")) %>%
    # Remove extra semicolons and whitespace
    mutate(Modifications = str_replace_all(Modifications, "\\s*;\\s*", ";")) %>%
    mutate(Modifications = str_replace_all(Modifications, ";{2,}", ";")) %>%
    mutate(Modifications = str_remove_all(Modifications, "^;|;$"))

  # Flag rows as oxidized or unoxidized
  raw_data <- raw_data %>%
    mutate(MOD = ifelse(Modifications == "" | is.na(Modifications) | Modifications == "NA",
                        "Unoxidized", "Oxidized"))

  # Extract the first modified residue and its position (if multiple mods, just get the first)
  raw_data <- raw_data %>%
    mutate(first_mod = str_extract(Modifications, "[A-Z](?=\\d+\\()"),  # extract residue
           first_pos = str_extract(Modifications, "(?<=\\D)\\d+(?=\\()"))  # extract position

  # Assign to new columns
  raw_data <- raw_data %>%
    mutate(ModPositionL = first_mod,
           ModPositionN = as.numeric(first_pos))

  return(raw_data)
}

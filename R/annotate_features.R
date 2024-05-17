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

  raw_data <- raw_data %>%
    rename(Mods = Modifications)

  raw_data <- raw_data %>%
    mutate(Modifications = Mods) %>%
    mutate(Modifications = gsub(";[A-Z]\\d+\\(Carbamidomethyl\\)", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\);", ";", Modifications)) %>%
    mutate(Modifications = gsub("[A-Z]\\d+\\(Carbamidomethyl\\)", "", Modifications))


  raw_data <- raw_data %>%
    mutate(Modifications = gsub("\\s*;\\s*", ";", Modifications)) %>%  # Remove spaces around semicolons
    mutate(Modifications = gsub(";{2,}", ";", Modifications)) %>%  # Replace multiple semicolons with a single one
    mutate(Modifications = gsub("^;|;$", "", Modifications))  # Remove leading and trailing semicolons

  raw_data <- raw_data %>%
    mutate(MOD = ifelse(is.na(Modifications) | Modifications == "", "Unoxidized",
                        ifelse(grepl("", Modifications), "Oxidized", "Unoxidized")))

  raw_data <- raw_data %>%
    mutate(ModPositionL = sub("^([A-Z]*).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = gsub(".*?([0-9]+).*", "\\1", Modifications)) %>%
    mutate(ModPositionN = ifelse(ModPositionN == Modifications, NA, ModPositionN)) %>%
    mutate(ModPositionN = as.numeric(ModPositionN)) %>%  # Ensure ModPositionN is numeric
    mutate(ModPosition = ifelse(is.na(ModPositionL) | ModPositionL == "", NA,
                                paste(ModPositionL, ModPositionN, sep = "")))

  return(raw_data)
}



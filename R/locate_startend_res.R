#' Locate the Start and End Residue of a Peptide (Step 8)
#'
#' @param raw_data Annotated data frame to map peptide locations on
#' protein sequences
#'
#' @return A data frame with with the start and ending peptide location
#' @export
#'
#' @examples new_df <- locate_startend_res(raw_data)
#' @aliases locate_startend_res
locate_startend_res <- function(raw_data, FASTA){

  # Merge raw_data with FASTA by UniprotID
  raw_data <- merge(raw_data, FASTA, by = "UniprotID", all.x = TRUE)

  # Locate the start and end of the sequence within the protein sequence
  index <- str_locate(raw_data$protein_sequence, raw_data$Sequence)
  raw_data <- cbind(raw_data, index)

  # Create peptide column
  raw_data$peptide <- paste(raw_data$start, "-", raw_data$end, sep = "")

  # Count the number of modifications
  raw_data$mod_count <- str_count(raw_data$Modifications, "\\(.*?\\)")
  raw_data$mod_count <- ifelse(raw_data$MOD == "Unoxidized", 0, raw_data$mod_count)

  # Convert ModPositionN and start to numeric
  raw_data$ModPositionN <- as.numeric(raw_data$ModPositionN)
  raw_data$start <- as.numeric(raw_data$start)

  # Calculate modified residue positions
  raw_data$mod_res <- ifelse(!is.na(raw_data$ModPositionN) & raw_data$ModPositionN > 0, raw_data$start + raw_data$ModPositionN - 1, NA)

  # Create the Res column
  raw_data$Res <- paste(raw_data$ModPositionL, raw_data$mod_res, sep = "")

  return(raw_data)
}

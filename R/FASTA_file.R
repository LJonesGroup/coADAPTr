#' Import the Corresponding FASTA File
#'
#' @param FASTA_path FASTA file path that is user defined in the File Explorer Prompt
#' @return FASTA file saved as a data frame
#' @export
#'
#' @examples HumanFASTA <- FASTA_file(file.fasta)
#' @aliases FASTA_file
FASTA_file <- function() {
  cat("Please select the FASTA file associated with your data: \n")
  FASTA_path <- file.choose()
  FASTA <- read.fasta(FASTA_path)
  return(FASTA)
}



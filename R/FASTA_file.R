#' FASTA_file
#' @param FASTA_path FASTA file path
#' @return FASTA file saved as a data frame
#' @export
#'
#' @examples HumanFASTA <- FASTA_file(file.fasta)
#' @aliases FASTA_file
FASTA_file <- function() {
  FASTA_path <- file.choose()
  FASTA <- read.fasta(FASTA_path)
  return(FASTA)
}


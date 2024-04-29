#' parse fasta
#'
#' @param fasta_in
#'
#' @return a condensed dataframe with the protein sequence and the uniprot ID
#' @export
#'
#' @examples filtered_FASTA <- parse_fasta(FASTA)
#' @aliases parse
parse_fasta <- function(fasta_in){
  fasta_in <- fasta_in %>%
    rename(protein_sequence = seq.text)
  fasta_in<- fasta_in %>%
    rename(MasterProteinAccessions = seq.name)

  fasta_in$UniprotID <- gsub("^.+\\|(\\w+)\\|.*$", "\\1", fasta_in$MasterProteinAccessions)

  return(fasta_in)

}

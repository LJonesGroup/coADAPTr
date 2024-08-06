
#' Prepare the Input TMT Data for Calculating the Extent of Modification by
#' Identifying the Sample vs Control Files and  Start and End of the Peptide Sequence
#'
#' @return Data Frames Containing Annotated TMT Data
#' @export
#'
#' @examples TMT_Prep()
#' @aliases TMT_Prep
TMT_Prep <- function() {
  #Identify Sample vs Control Files
  transformed_data <<- SampleControl(transformed_data)
  #Annotate Features in the Data
  annotated_data <<- annotate_features(transformed_data)
  #Parse the FASTA file for later manipulations
  FASTA<<- parse_fasta(FASTA)
  #Locate the Start and End Residues for Each Peptide
  mod_data_fasta_merged<<- locate_startend_res(annotated_data, FASTA)
}


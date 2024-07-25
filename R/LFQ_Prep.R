
#' Prepare the Input Data for Calculating the Extent of Modification by
#' Selecting and Renaming the Relevant Columns and IDentifying the Start and End
#' of the Peptide Sequence
#'
#' @return Data Frames Containing Annotated Data
#' @export
#'
#' @examples LFQ_Prep()
#' @aliases LFQ_Prep
LFQ_Prep <- function() {
  #Select Relevant Data
  selected_data <<- column_selectionLFQ(raw_data)
  #Rename Certain Columns and Split the Spectrum Files
  modified_data <<- rename_and_split_spectrum_files(selected_data)
  #Identify Sample vs Control Files
  modified_data <<- SampleControl(modified_data)
  #Annotate Features in the Data
  modified_data_annotated <<- annotate_features(modified_data)
  #Parse the FASTA file for later manipulations
  FASTA<<- parse_fasta(FASTA)
  #Locate the Start and End Residues for Each Peptide
  mod_data_fasta_merged<<- locate_startend_res(modified_data_annotated, FASTA)
}


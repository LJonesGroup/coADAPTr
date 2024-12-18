
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
  # Detect raw_data type based on finding a specific column name and choose
  # for the users

  columns_rawdata <- colnames(raw_data)

  # TRUE (Correct)
  #print("Q Value" %in% columns_rawdata)
  # False (So the point get converted to a space)
  #print("Q.Value" %in% columns_rawdata)

  if ("Delta Mass" %in% columns_rawdata & "Assigned Modifications" %in% columns_rawdata){
    # FragPipe DDA
    print("Extracting columns from psm.tsv file")
    selected_data <<- column_selectionLFQ(raw_data, "psm")
    #Rename Certain Columns and Split the Spectrum Files
    modified_data <<- rename_and_split_spectrum_files(selected_data, TRUE)
    #Identify Sample vs Control Files
    modified_data <<- SampleControl(modified_data)
    #Annotate Features in the Data
    modified_data_annotated <<- annotate_features(modified_data)
    #Parse the FASTA file for later manipulations
    FASTA<<- parse_fasta(FASTA)
    #Locate the Start and End Residues for Each Peptide
    mod_data_fasta_merged<<- locate_startend_res(modified_data_annotated, FASTA)

  }else if("Q Value" %in% columns_rawdata & "Modified Sequence" %in% columns_rawdata){
    # FragPipe DIA
    print("Extracting columns from report.tsv file")
    #Select Relevant Data
    selected_data <<- column_selectionLFQ(raw_data, "report")
    #Rename Certain Columns and Split the Spectrum Files
    modified_data <<- rename_and_split_spectrum_files(selected_data, TRUE)
    #Identify Sample vs Control Files
    modified_data <<- SampleControl(modified_data)
    #Annotate Features in the Data
    modified_data_annotated <<- annotate_features(modified_data)
    #Parse the FASTA file for later manipulations
    FASTA<<- parse_fasta(FASTA)
    #Locate the Start and End Residues for Each Peptide
    mod_data_fasta_merged<<- locate_startend_res(modified_data_annotated, FASTA)


  }else{
    #Proteome Discover
    #Select Relevant Data
    selected_data <<- column_selectionLFQ(raw_data, "PD")
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


}


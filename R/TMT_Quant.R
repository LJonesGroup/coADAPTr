
#' Prepare the Input Data TMT for Annotation and Quantification
#'
#' @return Data Frames Containing Quantified, Renamed, and Transformed TMT Columns
#' @export
#'
#' @examples TMT_Quant()
#' @aliases TMT_Quant
TMT_Quant <- function() {
  #Summarize and Calculate the Normalized TMT Abundances
  TMT_data <<- sum_and_calculate_abundancesTMT(raw_data)
  #Select the required Columns
  selected_data <<- column_selectionTMT(TMT_data)
  #Rename the TMT columns to reflect the conditions and sample types
  renamed_data <<- rename_columns_interactively(selected_data)
  #Transform the data to a format that can be used for further analysis
  transformed_data <<- transform_data(renamed_data)
  #Identify Sample vs Control Files
  transformed_data <<- SampleControl(transformed_data)
}


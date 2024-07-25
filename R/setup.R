

#' Prepare R by Installing the Necessary Packages; Eliminating Conflicts, Importing Data, and
#' Selecting an Output Folder
#'
#' @return The imported FASTA file, raw data, and the selected output folder
#' @export
#'
#' @examples setup()
#' @aliases setup
setup <- function() {

  before_beginning()

  raw_data <<- import_data()

  file_output <<- output_folder()

  FASTA <<- FASTA_file()
}


